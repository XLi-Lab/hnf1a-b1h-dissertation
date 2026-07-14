# HNF1A B1H binding-specificity analysis

Analysis code for an MSc dissertation characterising HNF1A DNA-binding specificity by bacterial one-hybrid (B1H) selection, deep mutational scanning and structural modelling.

## Scripts

- `hnf1a_b1h_pipeline.py`: Phase 1 pipeline. Raw Nanopore amplicon FASTQs through demultiplexing, motif mapping, enrichment scoring, tier assignment, and gene and tissue annotation.
- `analyse_af3_fig7.py`: AlphaFold 3 analysis. Per-residue pLDDT, the strongest predicted DNA-contact residues, model confidence metrics, and the 21-site confidence-versus-binding comparison.
- `analyse_cici_dms.py`: deep mutational scanning analysis. Fitness-scale validation, DNA-side and protein-side single-site maps, per-residue sensitivity, protein-by-DNA coupling, and coupling at the DNA-binding interface.

Figures are produced by separate plotting scripts, which are not included here. Every statistic reported in the dissertation is computed by the scripts above.

## Requirements

Python 3.10+. Install dependencies with:

```
pip install -r requirements.txt
```

## Running

Place the input files listed at the top of each script alongside it, then run:

```
python hnf1a_b1h_pipeline.py
python analyse_af3_fig7.py
python analyse_cici_dms.py
```

## Data availability

Input datasets are not included in this repository.

The deep mutational scanning datasets were generated in the Lehner laboratory using the TF-MAPS system and are unpublished. They are available on request, subject to permission.

The Phase 1 sequencing data, the AlphaFold 3 model outputs and the processed enrichment tables are available from the author on request.

The Phase 2 analysis pipeline will be added to this repository once sequencing is complete.

## Author

Reem Sayed Ahmed, King's College London, 2026. Supervised by Dr Xianghua Li.
