# HNF1A B1H binding-specificity analysis

Analysis code for an MSc dissertation characterising HNF1A DNA-binding specificity by bacterial one-hybrid (B1H) selection, deep mutational scanning and structural modelling.

Reem Sayed Ahmed, King's College London, 2026. Supervised by Dr Xianghua Li, Li Laboratory, Guy's & St Thomas'.

## Contents

- `hnf1a_b1h_pipeline.py` — Phase 1 pipeline: takes the raw Nanopore amplicon reads through demultiplexing, motif mapping and enrichment scoring to the tissue-annotated hit table.
- `analyse_cici_dms.py` — deep mutational scanning analysis: fitness-scale validation, the DNA-side and protein-side single-site maps, and the protein-by-DNA coupling analysis.
- `analyse_af3_fig7.py` — AlphaFold 3 analysis: per-residue confidence, predicted DNA contacts, and interface confidence against experimental binding across 21 sites.
- `folds_2026_06_22_14_59.zip` — AlphaFold 3 model outputs for the 21 modelled sites (seven per class).
- `requirements.txt` — Python dependencies.

## Setup

Python 3.10 or newer. Install dependencies:

    pip install -r requirements.txt

## Running the AlphaFold analysis

`analyse_af3_fig7.py` expects the model outputs at `af3/folds_2026_06_22_14_59/`. Unzip the archive into an `af3/` folder first:

    mkdir -p af3 && unzip folds_2026_06_22_14_59.zip -d af3/

## Input data

Some inputs are not held in this repository and are available from the author on request:

- Phase 1 raw reads (`Lib-only_pool_Tube1_raw_fastq.gz`, `DBD-Lib_pool_Tube2_raw_fastq.gz`), the designed motif list (`hnf1a_63_unique_motif.txt`), the BLAT hits (`blat_best_hits.json`) and the UCSC gene-overlap table (`BLAT_gene_overlaps.tsv`). `hnf1a_b1h_pipeline.py` reads these from an `./input/` folder.
- The processed Phase 1 enrichment table (`HNF1A_B1H_all_motifs_corrected_2.csv`), needed by the DMS and AlphaFold analysis scripts.
- The deep mutational scanning datasets were generated in the Lehner laboratory using the TF-MAPS system and are unpublished. They are available on request, subject to permission.

The figure-generation scripts and the Phase 2 pipeline (pending sequencing) are also available from the author.
