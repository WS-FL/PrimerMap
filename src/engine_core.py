#!/usr/bin/env python3
"""Design local genomic PCR/Sanger primers around supplied 20-nt guides.

Input coordinates are zero-based, half-open on the supplied plus-strand sequence.
All specificity checks are restricted to that complete supplied sequence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import primer3
import primer3.bindings
from Bio import __version__ as biopython_version
from Bio.Seq import Seq


GAP_MIN = 101
FLANK = 800
MAX_TH = 47.0
NUM_RETURN = 30
PRODUCT_MIN = 450
PRODUCT_OPT = 500
PRODUCT_MAX = 550


def configure(product_min: int = 450, product_opt: int = 500, product_max: int = 550) -> None:
    global PRODUCT_MIN, PRODUCT_OPT, PRODUCT_MAX
    if not (300 <= product_min <= product_opt <= product_max <= 1000):
        raise ValueError("PCR product range must satisfy 300 <= min <= optimum <= max <= 1000")
    PRODUCT_MIN, PRODUCT_OPT, PRODUCT_MAX = product_min, product_opt, product_max
CHEMISTRY = {
    "monovalent_salt_mM": 50.0,
    "divalent_Mg_mM": 1.5,
    "dNTP_mM": 0.6,
    "DNA_nM": 50.0,
    "tm_formula": "SantaLucia 1998 nearest-neighbor (Primer3 code 1)",
    "salt_corrections": "SantaLucia 1998 (Primer3 code 1)",
}
PARAMETERS = {
    "guide_length_nt": 20,
    "min_clear_gap_each_side_nt": GAP_MIN,
    "product_size_nt": {"min": 450, "opt": 500, "max": 550},
    "primer_size_nt": {"min": 18, "opt": 20, "max": 25},
    "primer_tm_C": {"min": 58.0, "opt": 60.0, "max": 62.0},
    "max_pair_tm_difference_C": 2.0,
    "primer_GC_percent": {"min": 40.0, "max": 60.0},
    "thermodynamic_thresholds_C": {
        "self_any": MAX_TH, "self_end": MAX_TH,
        "hairpin": MAX_TH, "pair_any": MAX_TH, "pair_end": MAX_TH,
    },
    "local_flank_each_side_nt": FLANK,
    "primer3_candidates_requested": NUM_RETURN,
    "specificity_scope": "perfect-match sites on both strands of supplied full genomic sequence only; no genome-wide assessment",
    "chemistry": CHEMISTRY,
}


def revcomp(s: str) -> str:
    return str(Seq(s).reverse_complement())


def matches(s: str, motif: str) -> list[int]:
    positions = []
    start = 0
    while True:
        p = s.find(motif, start)
        if p < 0:
            return positions
        positions.append(p)
        start = p + 1


def match_sites_both_strands(genome: str, oligo: str) -> list[dict]:
    sites = [{"strand": "+", "start0": p} for p in matches(genome, oligo)]
    sites.extend({"strand": "-", "start0": p} for p in matches(genome, revcomp(oligo)))
    return sites


def inward_facing_products(fsites: list[dict], rsites: list[dict], flen: int, rlen: int) -> list[dict]:
    products = []
    for f in fsites:
        for r in rsites:
            if f["strand"] == "+" and r["strand"] == "-" and f["start0"] < r["start0"]:
                products.append({"start0": f["start0"], "end0": r["start0"] + rlen,
                                 "length_nt": r["start0"] + rlen - f["start0"], "orientation": "F+/R-"})
            if r["strand"] == "+" and f["strand"] == "-" and r["start0"] < f["start0"]:
                products.append({"start0": r["start0"], "end0": f["start0"] + flen,
                                 "length_nt": f["start0"] + flen - r["start0"], "orientation": "R+/F-"})
    return products


def gc_percent(s: str) -> float:
    return 100.0 * (s.count("G") + s.count("C")) / len(s)


def thermo_kwargs() -> dict:
    return {"mv_conc": 50.0, "dv_conc": 1.5, "dntp_conc": 0.6, "dna_conc": 50.0}


def calc_thermo(left: str, right: str) -> dict:
    kw = thermo_kwargs()
    lh = primer3.calc_hairpin(left, **kw)
    rh = primer3.calc_hairpin(right, **kw)
    ld = primer3.calc_homodimer(left, **kw)
    rd = primer3.calc_homodimer(right, **kw)
    hd = primer3.calc_heterodimer(left, right, **kw)
    return {
        "calculated_forward_hairpin_tm_C": lh.tm if lh.structure_found else None,
        "calculated_reverse_hairpin_tm_C": rh.tm if rh.structure_found else None,
        "calculated_forward_homodimer_tm_C": ld.tm if ld.structure_found else None,
        "calculated_reverse_homodimer_tm_C": rd.tm if rd.structure_found else None,
        "calculated_pair_heterodimer_tm_C": hd.tm if hd.structure_found else None,
    }


def primer3_args(local: str, gs: int, ge: int) -> tuple[dict, dict]:
    target_start = gs - GAP_MIN
    target_end = ge + GAP_MIN
    if target_start < 0 or target_end > len(local):
        raise ValueError("Guide lies too close to a source boundary for >100-nt gaps")
    seq_args = {
        "SEQUENCE_ID": "local_guide_amplicon",
        "SEQUENCE_TEMPLATE": local,
        "SEQUENCE_TARGET": [target_start, target_end - target_start],
        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [
            [0, target_start, target_end, len(local) - target_end]
        ],
    }
    global_args = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_NUM_RETURN": NUM_RETURN,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MAX_SIZE": 25,
        "PRIMER_MIN_TM": 58.0,
        "PRIMER_OPT_TM": 60.0,
        "PRIMER_MAX_TM": 62.0,
        "PRIMER_MAX_DIFF_TM": 2.0,
        "PRIMER_MIN_GC": 40.0,
        "PRIMER_MAX_GC": 60.0,
        "PRIMER_PRODUCT_SIZE_RANGE": [[PRODUCT_MIN, PRODUCT_MAX]],
        "PRIMER_PRODUCT_OPT_SIZE": PRODUCT_OPT,
        "PRIMER_SALT_MONOVALENT": 50.0,
        "PRIMER_SALT_DIVALENT": 1.5,
        "PRIMER_DNTP_CONC": 0.6,
        "PRIMER_DNA_CONC": 50.0,
        "PRIMER_TM_FORMULA": 1,
        "PRIMER_SALT_CORRECTIONS": 1,
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
        "PRIMER_MAX_SELF_ANY_TH": MAX_TH,
        "PRIMER_MAX_SELF_END_TH": MAX_TH,
        "PRIMER_MAX_HAIRPIN_TH": MAX_TH,
        "PRIMER_PAIR_MAX_COMPL_ANY_TH": MAX_TH,
        "PRIMER_PAIR_MAX_COMPL_END_TH": MAX_TH,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_MAX_POLY_X": 5,
        "PRIMER_EXPLAIN_FLAG": 1,
    }
    return seq_args, global_args


def primer3_explain(result: dict) -> dict:
    return {k: v for k, v in result.items() if "EXPLAIN" in k or k == "PRIMER_ERROR"}


def parse_candidate(result: dict, i: int, full: str, offset: int, gs: int, ge: int) -> dict:
    flocal, flen = result[f"PRIMER_LEFT_{i}"]
    rlast, rlen = result[f"PRIMER_RIGHT_{i}"]
    fs, fe = offset + flocal, offset + flocal + flen
    rs, re = offset + rlast - rlen + 1, offset + rlast + 1
    forward = result[f"PRIMER_LEFT_{i}_SEQUENCE"].upper()
    reverse = result[f"PRIMER_RIGHT_{i}_SEQUENCE"].upper()
    ftm = float(result[f"PRIMER_LEFT_{i}_TM"])
    rtm = float(result[f"PRIMER_RIGHT_{i}_TM"])
    fgc, rgc = gc_percent(forward), gc_percent(reverse)
    fsites = match_sites_both_strands(full, forward)
    rsites = match_sites_both_strands(full, reverse)
    products = inward_facing_products(fsites, rsites, len(forward), len(reverse))
    thermo = calc_thermo(forward, reverse)
    calculated_ftm = primer3.calc_tm(forward, **thermo_kwargs(), tm_method="santalucia", salt_corrections_method="santalucia")
    calculated_rtm = primer3.calc_tm(reverse, **thermo_kwargs(), tm_method="santalucia", salt_corrections_method="santalucia")
    p3_thermo = {
        "forward_self_any_tm_C": result.get(f"PRIMER_LEFT_{i}_SELF_ANY_TH"),
        "forward_self_end_tm_C": result.get(f"PRIMER_LEFT_{i}_SELF_END_TH"),
        "forward_hairpin_tm_C": result.get(f"PRIMER_LEFT_{i}_HAIRPIN_TH"),
        "reverse_self_any_tm_C": result.get(f"PRIMER_RIGHT_{i}_SELF_ANY_TH"),
        "reverse_self_end_tm_C": result.get(f"PRIMER_RIGHT_{i}_SELF_END_TH"),
        "reverse_hairpin_tm_C": result.get(f"PRIMER_RIGHT_{i}_HAIRPIN_TH"),
        "pair_any_tm_C": result.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH"),
        "pair_end_tm_C": result.get(f"PRIMER_PAIR_{i}_COMPL_END_TH"),
    }
    candidate = {
        "primer3_candidate_index": i,
        "forward": {"sequence_5to3": forward, "binding_start0": fs, "binding_end0": fe,
                    "tm_C": ftm, "GC_percent": round(fgc, 2), "full_sequence_perfect_match_count_both_strands": len(fsites)},
        "reverse": {"sequence_5to3": reverse, "binding_start0": rs, "binding_end0": re,
                    "tm_C": rtm, "GC_percent": round(rgc, 2), "full_sequence_perfect_match_count_both_strands": len(rsites)},
        "product_start0": fs,
        "product_end0": re,
        "product_length_nt": re - fs,
        "left_clear_gap_nt": gs - fe,
        "right_clear_gap_nt": rs - ge,
        "tm_difference_C": round(abs(ftm - rtm), 3),
        "calculated_nn_tm_C": {"forward": calculated_ftm, "reverse": calculated_rtm},
        "primer3_pair_penalty": result.get(f"PRIMER_PAIR_{i}_PENALTY"),
        "thermodynamics_C": {**p3_thermo, **thermo},
        "full_sequence_match_sites": {"forward": fsites, "reverse": rsites},
        "full_sequence_inward_facing_perfect_match_product_count": len(products),
        "full_sequence_inward_facing_perfect_match_products": products[:5],
    }
    rejection = validate_candidate(candidate, full, gs, ge)
    candidate["accepted"] = not rejection
    candidate["rejection_reasons"] = rejection
    return candidate


def validate_candidate(c: dict, full: str, gs: int, ge: int) -> list[str]:
    reasons = []
    f, r = c["forward"], c["reverse"]
    fs, fe, rs, re = f["binding_start0"], f["binding_end0"], r["binding_start0"], r["binding_end0"]
    if not (0 <= fs < fe <= len(full) and 0 <= rs < re <= len(full)):
        return ["binding coordinates outside source sequence"]
    if full[fs:fe] != f["sequence_5to3"]:
        reasons.append("forward primer does not match plus strand at reported site")
    if full[rs:re] != revcomp(r["sequence_5to3"]):
        reasons.append("reverse primer does not reverse-complement reported plus-strand site")
    if fe - fs != len(f["sequence_5to3"]) or re - rs != len(r["sequence_5to3"]):
        reasons.append("binding interval length mismatch")
    if not (18 <= fe - fs <= 25 and 18 <= re - rs <= 25):
        reasons.append("primer size outside 18-25 nt")
    if gs - fe < GAP_MIN or rs - ge < GAP_MIN:
        reasons.append("guide-to-binding clear gap is not >100 nt on both sides")
    if not (PRODUCT_MIN <= re - fs <= PRODUCT_MAX):
        reasons.append(f"product size outside {PRODUCT_MIN}-{PRODUCT_MAX} nt")
    if not (58 <= f["tm_C"] <= 62 and 58 <= r["tm_C"] <= 62):
        reasons.append("primer Tm outside 58-62 C")
    if abs(f["tm_C"] - r["tm_C"]) > 2.0001:
        reasons.append("pair Tm difference >2 C")
    if abs(f["tm_C"] - c["calculated_nn_tm_C"]["forward"]) > 0.02 or abs(r["tm_C"] - c["calculated_nn_tm_C"]["reverse"]) > 0.02:
        reasons.append("Primer3 returned Tm disagrees with independent SantaLucia NN calculation")
    if not (40 <= f["GC_percent"] <= 60 and 40 <= r["GC_percent"] <= 60):
        reasons.append("primer GC outside 40-60%")
    if c["full_sequence_match_sites"]["forward"] != [{"strand": "+", "start0": fs}]:
        reasons.append("forward primer is not unique on both strands of supplied full sequence")
    if c["full_sequence_match_sites"]["reverse"] != [{"strand": "-", "start0": rs}]:
        reasons.append("reverse primer is not unique on both strands of supplied full sequence")
    if c["full_sequence_inward_facing_perfect_match_product_count"] != 1 or c["full_sequence_inward_facing_perfect_match_products"] != [
        {"start0": fs, "end0": re, "length_nt": re - fs, "orientation": "F+/R-"}
    ]:
        reasons.append("pair does not yield exactly the intended inward-facing perfect-match product in supplied full sequence")
    for key, value in c["thermodynamics_C"].items():
        if value is not None and value > MAX_TH:
            reasons.append(f"{key} > {MAX_TH} C")
    return reasons


def design_one(rec: dict) -> dict:
    gene = str(rec.get("gene", "unknown"))
    out = {"gene": gene, "source_file": rec.get("source_file"), "status": "failed",
           "selected_pairs": [], "candidates": [], "failure_reasons": []}
    try:
        full = re.sub(r"\s+", "", rec["sequence"]).upper()
        if not full or re.search(r"[^ACGTN]", full):
            raise ValueError("Source sequence must contain only A/C/G/T/N")
        gs, ge = int(rec["guide_start0"]), int(rec["guide_end0"])
        guide = str(rec["guide_sequence"]).upper()
        strand_raw = str(rec.get("guide_strand", "")).strip()
        strand = "+" if strand_raw in ("+", "+1", "1") else "-" if strand_raw in ("-", "-1") else strand_raw
        if not (0 <= gs < ge <= len(full) and ge - gs == 20 and len(guide) == 20):
            raise ValueError("Guide interval and sequence must both be exactly 20 nt within source")
        observed = full[gs:ge]
        if strand == "+" and observed != guide:
            raise ValueError("Plus-strand guide sequence does not match supplied interval")
        if strand == "-" and revcomp(observed) != guide:
            raise ValueError("Minus-strand guide sequence does not match supplied interval")
        if strand not in ("+", "-") and guide not in (observed, revcomp(observed)):
            raise ValueError("Guide sequence does not match either strand of supplied interval")
        if gs < GAP_MIN or len(full) - ge < GAP_MIN:
            raise ValueError("Guide too close to full-sequence boundary for required clear gaps")
        offset = max(0, gs - FLANK)
        end = min(len(full), ge + FLANK)
        local = full[offset:end]
        seq_args, glob_args = primer3_args(local, gs - offset, ge - offset)
        result = primer3.bindings.design_primers(seq_args, glob_args)
        out.update({"full_source_length_nt": len(full), "guide_start0": gs, "guide_end0": ge,
                    "guide_strand": strand, "guide_sequence_5to3": guide,
                    "local_template_start0": offset, "local_template_end0": end,
                    "primer3_global_args": glob_args,
                    "primer3_sequence_constraints_local": {
                        "SEQUENCE_TARGET": seq_args["SEQUENCE_TARGET"],
                        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": seq_args["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"],
                    },
                    "primer3_explain": primer3_explain(result),
                    "primer3_returned_pair_count": result.get("PRIMER_PAIR_NUM_RETURNED", 0)})
        if result.get("PRIMER_ERROR"):
            raise ValueError(f"Primer3: {result['PRIMER_ERROR']}")
        for i in range(result.get("PRIMER_PAIR_NUM_RETURNED", 0)):
            out["candidates"].append(parse_candidate(result, i, full, offset, gs, ge))
        valid = [c for c in out["candidates"] if c["accepted"]]
        valid.sort(key=lambda c: (abs(c["product_length_nt"] - PRODUCT_OPT), c["primer3_pair_penalty"]))
        if valid:
            first = valid[0]
            remainder = [c for c in valid[1:] if
                         (c["forward"]["sequence_5to3"], c["reverse"]["sequence_5to3"]) !=
                         (first["forward"]["sequence_5to3"], first["reverse"]["sequence_5to3"])]
            remainder.sort(key=lambda c: (
                -(int(c["forward"]["sequence_5to3"] != first["forward"]["sequence_5to3"]) +
                  int(c["reverse"]["sequence_5to3"] != first["reverse"]["sequence_5to3"])),
                abs(c["product_length_nt"] - PRODUCT_OPT), c["primer3_pair_penalty"]))
            chosen = [first] + remainder[:1]
            for label, c in zip(("primary", "backup"), chosen):
                assert not validate_candidate(c, full, gs, ge)
                c["selection"] = label
                out["selected_pairs"].append(c)
            out["status"] = "success"
        else:
            out["failure_reasons"].append("No candidate passed independent constraints and full-source two-strand uniqueness checks")
            if out["candidates"]:
                reasons = sorted({reason for c in out["candidates"] for reason in c["rejection_reasons"]})
                out["failure_reasons"].extend(reasons)
            else:
                out["failure_reasons"].append("Primer3 returned zero candidate pairs")
    except Exception as exc:
        out["failure_reasons"].append(f"{type(exc).__name__}: {exc}")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path(__file__).with_name("inputs.json"))
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("design_results.json"))
    args = parser.parse_args()
    raw = args.input.read_bytes()
    payload = json.loads(raw)
    records = payload["records"]
    if not isinstance(records, list):
        raise ValueError("inputs.json records must be a list")
    result = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256": hashlib.sha256(raw).hexdigest(),
        "software": {"primer3_py": primer3.__version__, "biopython": biopython_version},
        "parameters": PARAMETERS,
        "record_count": len(records),
        "results": [design_one(rec) for rec in records],
    }
    result["success_count"] = sum(r["status"] == "success" for r in result["results"])
    result["failed_count"] = len(records) - result["success_count"]
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(f"Designed {result['success_count']}/{len(records)} records; output: {args.output}")


if __name__ == "__main__":
    main()
