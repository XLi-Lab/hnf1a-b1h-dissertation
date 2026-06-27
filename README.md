# HNF1A B1H binding-specificity analysis

Analysis code for an MSc dissertation characterising HNF1A DNA-binding
specificity by bacterial one-hybrid (B1H) selection.

## Scripts

- `hnf1a_b1h_pipeline.py`: Phase 1 pipeline: raw Nanopore amplicon FASTQs
  through to enrichment scores, tier assignment, gene and tissue annotation,
  and the 3-panel figure.
- `analyse_af3_fig7.py`: reproduces the AlphaFold 3 numbers (Figure 7):
  per-residue pLDDT, the strongest DNA-contact residues, confidence metrics,
  and the 21-site confidence-versus-binding comparison.
- `analyse_cici_dms.py`: deep mutational scanning analysis: fitness-scale
  validation, single-site maps, per-residue sensitivity, and protein x DNA
  coupling.

## Requirements

Python 3.10+. Install dependencies with:

    pip install -r requirements.txt

## Running

Place the input files in the locations noted at the top of each script, then run:

    python hnf1a_b1h_pipeline.py
    python analyse_af3_fig7.py
    python analyse_cici_dms.py

## Data availability

The input datasets are not included in this repository. The deep mutational
scanning datasets are unpublished and available from the Li Lab on request.

## Author

Reem Sayed Ahmed, King's College London, 2026.
