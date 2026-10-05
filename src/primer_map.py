#!/usr/bin/env python3
"""Local Primer3/Sanger mapping from sg-marked GenBank files. No web server."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import warnings
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from reportlab.lib import colors
from reportlab.lib.pagesizes import A3, landscape
from reportlab.pdfgen import canvas

import engine_core as core


def resource_path(name: str) -> Path:
    root = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent.parent))
    return root / "reference_cache" / name


def read_gb(path: Path) -> SeqRecord:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rec = SeqIO.read(path, "genbank")
    seq = str(rec.seq).upper()
    if not seq or set(seq) - set("ACGTN"):
        raise ValueError("GenBank DNA must contain only A/C/G/T/N")
    rec.seq = Seq(seq)
    return rec


def q(feature: SeqFeature, name: str) -> str:
    return " ".join(feature.qualifiers.get(name, []))


def gene_name(path: Path, rec: SeqRecord) -> str:
    first = re.match(r"[A-Za-z][A-Za-z0-9-]*", path.stem)
    if first:
        return first.group(0).upper()
    for f in rec.features:
        if q(f, "gene"):
            return q(f, "gene").upper()
    return re.sub(r"\W+", "_", rec.name).upper()


def guides(rec: SeqRecord) -> tuple[list[dict], list[str]]:
    seq = str(rec.seq).upper()
    found, ignored = [], []
    for f in rec.features:
        if f.type not in ("primer", "primer_bind", "misc_feature"):
            continue
        label, note = q(f, "label"), q(f, "note")
        if re.search(r"sequencing|(?:^|\W)r?ww\d+", label, re.I):
            continue
        m = re.search(r"(?:^|\b)sequence\s*:\s*([ACGTacgt]{20})(?![ACGTacgt])", note)
        if not m:
            if f.type == "misc_feature" and re.search(r"\bsg(?:RNA)?\s*\d*\b", label, re.I) and len(f.location) == 20:
                strand = f.location.strand or 1
                s = seq[int(f.location.start):int(f.location.end)]
                mseq = s if strand == 1 else core.revcomp(s)
            else:
                continue
        else:
            mseq = m.group(1).upper()
        a, b = int(f.location.start), int(f.location.end)
        strand = f.location.strand or 1
        if b - a != 20:
            ignored.append(f"{label or f.type}: annotation length is not 20 nt")
            continue
        observed = seq[a:b] if strand == 1 else core.revcomp(seq[a:b])
        if observed != mseq:
            ignored.append(f"{label or f.type}: marked oligo does not exactly match its annotated strand/site")
            continue
        found.append({"start0": a, "end0": b, "strand": strand, "sequence": mseq, "label": label})
    unique = {(g["start0"], g["end0"], g["strand"], g["sequence"]): g for g in found}
    return sorted(unique.values(), key=lambda g: g["start0"]), ignored


def source_gene_id(rec: SeqRecord) -> str | None:
    for f in rec.features:
        if f.type == "gene":
            m = re.search(r"GeneID\s*:\s*(\d+)", q(f, "note") + " " + q(f, "db_xref"))
            if m:
                return m.group(1)
    return None


def ncbi_text(url: str, max_bytes: int = 8_000_000) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "PrimerMap/1.0 (local desktop tool)"})
    with urllib.request.urlopen(req, timeout=35) as response:
        data = response.read(max_bytes + 1)
    if len(data) > max_bytes:
        raise ValueError("NCBI reference response too large; provide a shorter annotated GB")
    return data.decode("utf-8")


def ncbi_reference(rec: SeqRecord, gene: str, out_cache: Path) -> tuple[SeqRecord, str]:
    # An included, previously independently verified GRCh38 fragment makes the
    # known MFSD12 boundary case work even without network access.
    cached = resource_path("MFSD12_extended_annotated.gb")
    if gene == "MFSD12" and cached.exists():
        ref = read_gb(cached)
        if str(rec.seq) in str(ref.seq):
            return ref, "NCBI GRCh38 NC_000019.10; bundled, sequence-verified MFSD12 reference"
    gene_id = source_gene_id(rec)
    if not gene_id:
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(
            {"db": "gene", "term": f"{gene}[Gene Name] AND Homo sapiens[Organism]", "retmode": "json"})
        ids = json.loads(ncbi_text(search_url))["esearchresult"]["idlist"]
        if len(ids) != 1:
            raise ValueError(f"Cannot resolve {gene} uniquely to NCBI GeneID; matches: {ids[:5]}")
        gene_id = ids[0]
    es_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(
        {"db": "gene", "id": gene_id, "retmode": "json"})
    data = json.loads(ncbi_text(es_url))
    summary = data["result"][gene_id]
    if summary.get("name", "").upper() != gene.upper():
        aliases = [x.upper() for x in summary.get("otheraliases", "").split(", ")]
        if gene.upper() not in aliases:
            raise ValueError(f"NCBI GeneID {gene_id} resolves to {summary.get('name')}, not {gene}")
    infos = [x for x in summary.get("genomicinfo", []) if x.get("chraccver", "").startswith("NC_")]
    if not infos:
        raise ValueError("NCBI gene summary has no RefSeq chromosome accession")
    info = infos[0]
    accession = info["chraccver"]
    start = max(1, min(info["chrstart"], info["chrstop"]) + 1 - 1200)
    end = max(info["chrstart"], info["chrstop"]) + 1 + 1200
    if end - start + 1 > 1_000_000:
        raise ValueError("Reference interval exceeds 1 Mb; provide a shorter annotated GB")
    cache_path = out_cache / f"{accession}_{start}_{end}.gb"
    if cache_path.exists():
        ref = read_gb(cache_path)
    else:
        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(
            {"db": "nuccore", "id": accession, "seq_start": start, "seq_stop": end,
             "rettype": "gbwithparts", "retmode": "text"})
        raw = ncbi_text(fetch_url)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ref = SeqIO.read(io.StringIO(raw), "genbank")
        if len(ref.seq) != end - start + 1:
            raise ValueError("NCBI returned a reference length different from the requested coordinates")
        cache_path.write_text(raw)
        time.sleep(0.35)
    return ref, f"NCBI {accession}:{start}-{end}; GeneID {gene_id}"


def feature_exons(rec: SeqRecord) -> list[SeqFeature]:
    return [f for f in rec.features if f.type == "exon"]


def prepare_record(source: SeqRecord, guide: dict, gene: str, allow_network: bool, cache: Path,
                   force_reference: bool = False) -> tuple[SeqRecord, dict, str]:
    old = str(source.seq).upper()
    enough_flank = min(guide["start0"], len(old) - guide["end0"]) >= 101
    if feature_exons(source) and enough_flank and not force_reference:
        return copy.deepcopy(source), guide.copy(), "source GB exon annotations"
    if not allow_network and gene != "MFSD12":
        raise ValueError("Exon annotation or sufficient flanks missing; enable NCBI retrieval")
    ref, provenance = ncbi_reference(source, gene, cache)
    plus = str(ref.seq).upper()
    pos_plus = plus.find(old)
    pos_minus = core.revcomp(plus).find(old)
    if (pos_plus >= 0) == (pos_minus >= 0):
        raise ValueError("Full GB sequence is absent or ambiguous in reference; no safe coordinate mapping")
    oriented_ref = ref if pos_plus >= 0 else ref.reverse_complement(id=True, name=True, description=True)
    oriented = str(oriented_ref.seq).upper()
    offset = oriented.find(old)
    if offset < 0 or oriented.find(old, offset + 1) >= 0:
        raise ValueError("Full GB sequence does not map uniquely to the reference")
    new = copy.deepcopy(source)
    new.seq = Seq(oriented)
    new.id = re.sub(r"\W+", "_", gene + "_extended")[:16]
    new.name = new.id
    new.description = f"{gene} reference-verified extended genomic DNA"
    new.features = [copy.deepcopy(f) for f in source.features]
    for f in new.features:
        f.location = f.location + offset
    mapped = guide.copy()
    mapped["start0"] += offset
    mapped["end0"] += offset
    if oriented[mapped["start0"]:mapped["end0"]] != (mapped["sequence"] if mapped["strand"] == 1 else core.revcomp(mapped["sequence"])):
        raise ValueError("Guide did not map to the verified reference")
    source_exons = feature_exons(source)
    if not source_exons:
        cached_exons = feature_exons(oriented_ref)
        if cached_exons:
            for exon in cached_exons:
                new.features.append(copy.deepcopy(exon))
            provenance += "; exon blocks from verified bundled reference"
            cached_exons = feature_exons(new)
        if cached_exons:
            new.annotations["comment"] = str(new.annotations.get("comment", "")) + f"\nPrimerMap verified complete source sequence at offset {offset+1}; {provenance}."
            return new, mapped, provenance
        enst = set(re.findall(r"ENST\d+\.\d+", str(source.annotations) + " " + " ".join(q(f, "note") for f in source.features)))
        mrnas = []
        for f in oriented_ref.features:
            if f.type != "mRNA" or len(f.location.parts) < 2:
                continue
            if not any(int(p.start) <= mapped["start0"] and int(p.end) >= mapped["end0"] for p in f.location.parts):
                continue
            tid = q(f, "transcript_id")
            if not tid.startswith("NM_"):
                continue
            xrefs = set(q(f, "db_xref").split())
            score = (bool(enst.intersection({x.replace("Ensembl:", "") for x in xrefs})), len(f.location.parts), tid)
            mrnas.append((score, f))
        if not mrnas:
            raise ValueError("Reference has no verified multi-exon NM transcript containing the whole guide")
        mrnas.sort(key=lambda x: x[0], reverse=True)
        best = mrnas[0][1]
        tid = q(best, "transcript_id")
        for i, part in enumerate(sorted(best.location.parts, key=lambda p: int(p.start)), 1):
            new.features.append(SeqFeature(FeatureLocation(int(part.start), int(part.end), strand=part.strand),
                type="exon", qualifiers={"label": [f"{gene} {tid} exon {i}"], "number": [str(i)],
                                         "transcript_id": [tid], "note": [provenance]}))
        provenance += f"; exon transcript {tid}"
    new.annotations["comment"] = str(new.annotations.get("comment", "")) + f"\nPrimerMap verified complete source sequence at offset {offset+1}; {provenance}."
    return new, mapped, provenance


def add_primer_features(rec: SeqRecord, g: dict, pair: dict) -> SeqRecord:
    out = copy.deepcopy(rec)
    f, r = pair["forward"], pair["reverse"]
    out.features.append(SeqFeature(FeatureLocation(g["start0"], g["end0"], strand=g["strand"]),
        type="misc_feature", qualifiers={"label": ["sgRNA 20 nt"], "note": [g["sequence"]]}))
    for name, p, strand in (("Sequencing F", f, 1), ("Sequencing R", r, -1)):
        out.features.append(SeqFeature(FeatureLocation(p["binding_start0"], p["binding_end0"], strand=strand),
            type="primer_bind", qualifiers={"label": [name], "note": [f"sequence: {p['sequence_5to3']}; Tm {p['tm_C']:.2f} C"]}))
    out.features.append(SeqFeature(FeatureLocation(pair["product_start0"], pair["product_end0"]),
        type="misc_feature", qualifiers={"label": [f"PCR product {pair['product_length_nt']} bp"]}))
    out.annotations["molecule_type"] = "DNA"
    return out


def write_pdf(path: Path, rows: list[dict], failures: list[dict]) -> None:
    W, H = landscape(A3)
    c = canvas.Canvas(str(path), pagesize=(W, H))
    c.setTitle("PrimerMap sequencing primer report")
    dark = colors.HexColor("#18334c")
    blue = colors.HexColor("#1871b5")
    purple = colors.HexColor("#8c46aa")
    red = colors.HexColor("#d64343")
    gray = colors.HexColor("#8396a9")
    def label(x, y, s, size=8, color=dark, bold=False):
        c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
        c.setFillColor(color)
        c.drawString(x, y, str(s)[:125])
    def arrow(a, b, y, color, forward=True):
        c.setStrokeColor(color); c.setLineWidth(3); c.line(a, y, b, y)
        x = b if forward else a
        d = -1 if forward else 1
        p = c.beginPath(); p.moveTo(x, y); p.lineTo(x+7*d, y+4); p.lineTo(x+7*d, y-4); p.close()
        c.setFillColor(color); c.drawPath(p, fill=1, stroke=0)
    label(35,H-40,"Sequencing primer map",21,dark,True)
    label(35,H-57,f"{len(rows)} designed guide(s)  |  {len(failures)} failed  |  arrows show primer extension",9,gray)
    y = H-90
    for row in rows:
        if y < 80:
            c.showPage(); y=H-45
        g,p = row["guide"],row["primary"]
        x0,x1=150,W-170
        v0=max(0,p["product_start0"]-75);v1=min(row["length"],p["product_end0"]+75)
        scale=lambda z:x0+(z-v0)/(v1-v0)*(x1-x0)
        label(35,y,row["gene"],10,dark,True)
        label(35,y-14,f"{p['product_length_nt']} bp",8,gray)
        c.setStrokeColor(dark);c.setLineWidth(1);c.line(x0,y,x1,y)
        for a,b in row["exons"]:
            if b<=v0 or a>=v1: continue
            c.setFillColor(colors.HexColor("#9caebe")); c.rect(scale(max(a,v0)),y-4,max(1,scale(min(b,v1))-scale(max(a,v0))),8,fill=1,stroke=0)
        f,r=p["forward"],p["reverse"]
        arrow(scale(f["binding_start0"]),scale(f["binding_end0"]),y+10,blue,True)
        arrow(scale(r["binding_start0"]),scale(r["binding_end0"]),y+10,purple,False)
        gx=scale((g["start0"]+g["end0"])/2)
        c.setFillColor(red);c.circle(gx,y,4,fill=1,stroke=0)
        label(W-155,y+2,f"gap {p['left_clear_gap_nt']}/{p['right_clear_gap_nt']}",8,dark)
        label(W-155,y-10,"F     sg     R",7,gray)
        y-=25
    label(35,38,"Exon = gray; forward = blue; sgRNA = red; reverse = purple. Each row uses its own amplicon scale.",8,gray)
    label(35,25,"Specificity: exact matches within supplied/extended gDNA only; no whole-genome screen.",8,gray)
    c.showPage()
    for row in rows:
        g,p=row["guide"],row["primary"]
        label(35,H-40,f"{row['gene']}  |  guide {row['guide_index']}",18,dark,True)
        label(35,H-60,f"Source: {row['source_name']}  |  {row['provenance']}",8,gray)
        label(35,H-88,f"sgRNA: {g['sequence']}  ({'+' if g['strand']==1 else '-'} strand), local bases {g['start0']+1}-{g['end0']}",10,red)
        label(35,H-110,f"Forward 5'-{p['forward']['sequence_5to3']}-3'   Tm {p['forward']['tm_C']:.2f} C",11,blue,True)
        label(35,H-132,f"Reverse 5'-{p['reverse']['sequence_5to3']}-3'   Tm {p['reverse']['tm_C']:.2f} C",11,purple,True)
        label(35,H-154,f"Amplicon {p['product_length_nt']} bp; guide-to-binding gaps {p['left_clear_gap_nt']} / {p['right_clear_gap_nt']} bp",10)
        label(35,H-176,f"Product bases {p['product_start0']+1}-{p['product_end0']}; F {p['forward']['binding_start0']+1}-{p['forward']['binding_end0']}; R {p['reverse']['binding_start0']+1}-{p['reverse']['binding_end0']}",9)
        label(35,H-198,f"Product range used {row['product_range'][0]}-{row['product_range'][2]} bp; guide within exon: {row['in_exon']}",9)
        if row.get("backup"):
            b=row["backup"]
            label(35,H-228,f"Backup: F {b['forward']['sequence_5to3']} / R {b['reverse']['sequence_5to3']} / {b['product_length_nt']} bp",9)
        label(35,60,"Primer3-py; SantaLucia NN; 50 mM monovalent, 1.5 mM Mg, 0.6 mM dNTP, 50 nM DNA.",8,gray)
        label(35,45,"No whole-genome specificity check. Verify template/genotype experimentally before ordering.",8,gray)
        c.showPage()
    if failures:
        label(35,H-40,"Files requiring attention",18,dark,True)
        y=H-75
        for x in failures:
            label(35,y,f"{x['source_name']}: {x['reason']}",8,red)
            y-=20
            if y<40: c.showPage();y=H-45
        c.showPage()
    c.save()


def process_file(path: Path, dest: Path, cache: Path, allow_network: bool, broad_fallback: bool) -> tuple[list[dict], list[dict]]:
    original = read_gb(path)
    gene = gene_name(path, original)
    marked, ignored = guides(original)
    if not marked:
        return [], [{"source_name": path.name, "reason": "No exact 20-nt sgRNA annotation found", "ignored": ignored}]
    successes, failures = [], []
    for index,g0 in enumerate(marked,1):
        try:
            rec,g,provenance=prepare_record(original,g0,gene,allow_network,cache)
            def run(working: SeqRecord, guide: dict, size: tuple[int,int,int]):
                core.configure(*size)
                return core.design_one({"gene":gene,"source_file":str(path),"sequence":str(working.seq),
                    "guide_start0":guide["start0"],"guide_end0":guide["end0"],
                    "guide_sequence":guide["sequence"],"guide_strand":guide["strand"]})
            ranges=[(450,500,550)] + ([(450,500,600),(350,500,700)] if broad_fallback else [])
            used=ranges[0]
            result=run(rec,g,used)
            for trial in ranges[1:]:
                if result["status"]=="success": break
                used=trial
                result=run(rec,g,used)
            if result["status"]!="success" and min(g["start0"],len(rec.seq)-g["end0"])<800 and allow_network:
                rec,g,provenance=prepare_record(original,g0,gene,allow_network,cache,force_reference=True)
                for trial in ranges:
                    used=trial
                    result=run(rec,g,used)
                    if result["status"]=="success": break
            if result["status"]!="success":
                raise ValueError("No pair passed constraints: " + "; ".join(result["failure_reasons"][:3]))
            exons=[(int(f.location.start),int(f.location.end)) for f in feature_exons(rec)]
            if not exons:
                raise ValueError("No verified exon features: position map would be incomplete")
            in_exon=any(a<=g["start0"] and b>=g["end0"] for a,b in exons)
            if not in_exon:
                raise ValueError("Guide is not fully inside any verified exon")
            p=result["selected_pairs"][0]
            b=result["selected_pairs"][1] if len(result["selected_pairs"])>1 else None
            annotated=add_primer_features(rec,g,p)
            source_tag=hashlib.sha256(path.read_bytes()).hexdigest()[:8]
            outfile=dest / f"{gene}_sg{index}_{source_tag}_sequencing_annotated.gb"
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                SeqIO.write(annotated,outfile,"genbank")
            reread=read_gb(outfile)
            if str(reread.seq)!=str(rec.seq) or not any(f.type=="exon" for f in reread.features):
                raise ValueError("Saved GenBank failed sequence/exon readback")
            successes.append({"gene":gene,"guide_index":index,"guide":g,"primary":p,"backup":b,
                "length":len(rec.seq),"exons":exons,"in_exon":in_exon,"product_range":used,
                "provenance":provenance,"source_name":path.name,"source_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
                "annotated_gb":outfile.name,"ignored_markers":ignored})
        except Exception as exc:
            failures.append({"source_name":path.name,"gene":gene,"guide_index":index,"reason":f"{type(exc).__name__}: {exc}","ignored":ignored})
    return successes,failures


def main() -> int:
    ap=argparse.ArgumentParser(description="Design local sequencing primers from sg-marked GenBank files")
    ap.add_argument("--input",type=Path,required=True,help=".gb/.gbk file or folder")
    ap.add_argument("--output",type=Path,required=True,help="output directory")
    ap.add_argument("--offline",action="store_true",help="never use NCBI to repair missing annotations/flanks")
    ap.add_argument("--strict-500",action="store_true",help="never widen 450-550 bp to 350-700 bp")
    args=ap.parse_args()
    paths=sorted(p for p in (args.input.iterdir() if args.input.is_dir() else [args.input])
                 if p.is_file() and p.suffix.lower() in (".gb",".gbk"))
    if not paths:
        ap.error("No .gb or .gbk files found")
    args.output.mkdir(parents=True,exist_ok=True)
    cache=args.output/"reference_cache";cache.mkdir(exist_ok=True)
    all_rows,all_fail=[],[]
    for i,path in enumerate(paths,1):
        print(f"[{i}/{len(paths)}] {path.name}",flush=True)
        try:
            rows,fail=process_file(path,args.output,cache,not args.offline,not args.strict_500)
        except Exception as exc:
            rows=[];fail=[{"source_name":path.name,"reason":f"{type(exc).__name__}: {exc}"}]
        all_rows+=rows;all_fail+=fail
        print(f"  designed {len(rows)}; attention {len(fail)}",flush=True)
    fields=["gene","guide","forward","reverse","f_tm_C","r_tm_C","amplicon_bp","left_gap_bp","right_gap_bp","product_range","annotated_gb","source_file"]
    with (args.output/"Primer_pairs.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields);w.writeheader()
        for row in all_rows:
            p=row["primary"]
            w.writerow({"gene":row["gene"],"guide":row["guide"]["sequence"],
                "forward":p["forward"]["sequence_5to3"],"reverse":p["reverse"]["sequence_5to3"],
                "f_tm_C":p["forward"]["tm_C"],"r_tm_C":p["reverse"]["tm_C"],
                "amplicon_bp":p["product_length_nt"],"left_gap_bp":p["left_clear_gap_nt"],
                "right_gap_bp":p["right_clear_gap_nt"],"product_range":f"{row['product_range'][0]}-{row['product_range'][2]}",
                "annotated_gb":row["annotated_gb"],"source_file":row["source_name"]})
    (args.output/"Review.json").write_text(json.dumps({"created_utc":datetime.now(timezone.utc).isoformat(),
        "input_count":len(paths),"success_count":len(all_rows),"failure_count":len(all_fail),
        "specificity_scope":"supplied full gDNA both strands only; no genomewide screen",
        "rows":all_rows,"failures":all_fail},ensure_ascii=False,indent=2)+"\n")
    with (args.output/"Needs_review.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=["source_name","gene","guide_index","reason"]);w.writeheader()
        for f in all_fail:w.writerow({k:f.get(k,"") for k in w.fieldnames})
    write_pdf(args.output/"Primer_map.pdf",all_rows,all_fail)
    print(f"Done: {len(all_rows)} guide(s), {len(all_fail)} attention item(s) -> {args.output}",flush=True)
    return 0 if not all_fail else 2


if __name__=="__main__":
    raise SystemExit(main())
