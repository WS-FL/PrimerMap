# GenBank input format

PrimerMap accepts one GenBank record per `.gb` or `.gbk` file. The record must contain genomic DNA (`A`, `C`, `G`, `T`, and optionally `N`) and at least one marked 20-nt sgRNA. The file name should start with the gene symbol, for example `SLC15A3 (ENSG00000110446).gb`; PrimerMap uses that first token as the output gene name.

## How PrimerMap finds a guide

**Recommended format for new files:** mark each guide as a `misc_feature` with a label beginning `sgRNA` or `sg` followed by an optional number, such as `sgRNA 1` or `sg2`. Its feature location must span exactly 20 bases. The GenBank location provides the position and strand; PrimerMap reads the 20-nt guide directly from the annotated sequence, reverse-complementing the interval when the feature is on the reverse strand. The PAM is **not** part of this 20-nt interval. A PAM or score may be added to the label or note for your own reference, but neither is used to locate the guide or design primers.

For example, a guide on bases 121–140 of the GenBank record (GenBank coordinates are 1-based and inclusive) can be represented as:

```text
     misc_feature    121..140
                     /label="sgRNA 1"
```

For a reverse-strand guide at the same bases, use `complement(121..140)`. In SnapGene or Benchling, create a 20-base feature named `sgRNA 1` and check the arrow direction and underlying bases before exporting a GenBank file. Each additional guide needs its own 20-base feature.

**Legacy compatibility:** PrimerMap also recognizes a 20-nt `primer` or `primer_bind` annotation when its `/note` contains `sequence: <20 A/C/G/T bases>`. It checks that sequence against the annotated interval in the indicated strand. This supports older files in which guides were saved as primers. Primer annotations whose labels include `sequencing` or a `ww#`/`rww#` index are excluded. Other 20-nt primer annotations with a matching `sequence:` note can still be interpreted as guides, so use the explicit `sgRNA` feature format for new inputs and remove unrelated primer annotations when using legacy files.

If multiple valid guide annotations exist, PrimerMap processes each distinct interval. Its output `sg1`, `sg2`, etc. indexes follow **position from left to right in the supplied record**, regardless of the number written in a feature label. Identical duplicate annotations at the same location, strand, and sequence are merged.

## Sequence context and checks

The genomic record needs enough sequence on **both** sides of every guide to place primer-binding sites more than 100 bases away from the entire guide. Longer flanks may be needed to find a pair satisfying the product-size and Tm constraints. Retain exon annotations if available: PrimerMap requires a verified exon containing the complete guide for its position map. If exons or flanks are missing, the app may retrieve a longer GRCh38/RefSeq record from NCBI and will transfer annotations only after exact, unique sequence matching. If it cannot verify that mapping, the target appears in `Needs_review.csv`.

The app does not infer a guide from a gene name, a PAM alone, a 20-base sequence in free text, or a bare unannotated GenBank sequence. It does not perform a genome-wide primer specificity check. It checks perfect-match primer binding on both strands of the supplied full sequence.

## What to verify before a batch run

1. Open each input in SnapGene or Benchling and confirm that the 20-base `sgRNA` feature covers the intended bases, excludes the PAM, and has the intended strand.
2. Confirm that the guide sits in the intended exon and that the file contains usable flanking genomic sequence on both sides.
3. Check `Primer_map.pdf` and `Review.json` after the run for the reported guide sequence, strand, local coordinates, exon, primer positions, and product length. Files without an accepted guide or primer pair are listed in `Needs_review.csv`.
