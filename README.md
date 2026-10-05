# PrimerMap for Mac

Developed by Boyuan Wang with implementation assistance from OpenAI Codex.

PrimerMap designs genomic PCR/Sanger primers around 20-nt sgRNA targets annotated in GenBank `.gb` or `.gbk` files. It runs locally without a browser, localhost server, or separate Python/Primer3 installation.

See [GenBank input format](INPUT_FORMAT.md) for the exact guide annotation, strand, and sequence-context requirements. For new files, mark each 20-nt guide as a `misc_feature` labeled `sgRNA 1`, `sgRNA 2`, and so on. The app uses the feature's location and strand to find the guide; it never guesses the target from the gene name alone.

Open `PrimerMap.app`, choose an input file or folder and an output folder, then click **Design Primers**. The results include:

- `Primer_pairs.csv`: forward/reverse primers, Tm, amplicon length, and guide-to-primer gaps.
- `*_sequencing_annotated.gb`: the original sequence and exon annotations with added sgRNA, primer, and PCR-product features.
- `Primer_map.pdf`: an overview map followed by target-level details.
- `Review.json` and `Needs_review.csv`: processing details and inputs requiring review.

The default product range is 450–550 bp. When no pair passes, PrimerMap tries 450–600 bp and then 350–700 bp. Check **Strict amplicon size** to keep the 450–550 bp range only. Both primer-binding sites must be more than 100 bp from the entire sgRNA; primer Tm must be 58–62 °C. Exact-match checks cover both strands of the supplied full genomic record. PrimerMap does not perform a genome-wide specificity screen.

When exons or flanking sequence are missing, PrimerMap retrieves a GRCh38/RefSeq record from NCBI. It requires the complete input sequence to match that record uniquely and exactly before moving annotations and restoring exon features. If it cannot verify the mapping, it puts the target in `Needs_review.csv` instead of guessing. A previously verified MFSD12 extension is bundled for offline use.

Inputs can mark a guide as a `primer`/`primer_bind` feature with `/note="sequence: <20 nt>"`, or as a `misc_feature` with an explicit `sgRNA` label and exact 20-nt span. Existing sequencing primers and `ww#`/`rww#` annotations are ignored as guide candidates. Distinct, valid guides in one input are processed separately.

This build targets Apple Silicon Macs. It is locally signed but not Apple-notarized. Other macOS versions and Gatekeeper behavior have not been tested.

Command-line use: `PrimerMap.app/Contents/Resources/engine/PrimerMapEngine --input /path/to/file.gb --output /path/to/results [--offline] [--strict-500]`.

## Build from source

On an Apple Silicon Mac with Python 3 and Xcode command-line tools:

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
zsh build_mac.sh
```

The resulting app and distributable ZIP are written to `dist/`. The build signs the app locally with an ad-hoc signature; it does not notarize it. Download a ready-to-use ZIP from the repository's Releases page if you do not need to build from source.
