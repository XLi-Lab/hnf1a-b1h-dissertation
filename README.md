# hnf1a-b1h-dissertation

Analysis code for the MSc dissertation *Dissecting the Determinants of HNF1A DNA-Binding
Specificity: A Bacterial One-Hybrid Screen of Native Promoters Integrated with Deep
Mutational Scanning of the POU-Homeodomain* (King's College London).

These are the analysis pipelines only. No plotting code is deposited. Each script is
standalone and prints every statistic quoted in the corresponding Results section,
including the assumption checks and the multiple-comparison corrections.

## Contents

| Script | Covers | Results section |
| --- | --- | --- |
| `hnf1a_b1h_pipeline.py` | Phase 1 bacterial one-hybrid: demultiplexing, counting, log2FC, tissue annotation, tier assignment | 3.1 to 3.6 |
| `analyse_consensus_distance.py` | Tier 1 consensus by majority vote and Hamming distance by functional class | 3.8 |
| `analyse_af3_confidence.py` | AlphaFold 3 model confidence for the 21 modelled sites | 3.10 |
| `analyse_cici_dms.py` | Deep mutational scan: site variants, coding variants, per-residue fitness, protein by DNA doubles | 3.11 to 3.15 |

## Statistical policy

The scripts follow the rule set out in Methods 2.18.

Comparisons of more than two groups run an omnibus test first. Shapiro-Wilk within each
group and Levene across groups decide the route: a one-way ANOVA with Tukey HSD where both
support it, otherwise Kruskal-Wallis. A post-hoc pairwise test is run only where the omnibus
is significant. Where a fixed number of planned two-group comparisons is made instead, the
p values are Bonferroni corrected for that number, and both the raw and the adjusted value
are printed so the correction is auditable.

There are four Bonferroni families, each of size two: silent versus missense and missense
versus nonsense; interface versus the rest and recognition helix versus all others; DNA core
versus flanks and protein interface versus elsewhere; and the AlphaFold three-class omnibus
with the pooled binder comparison. Single tests, such as the residue-shuffling permutation
test and the mutational-cost Spearman correlation, take no correction.

`analyse_consensus_distance.py` additionally repeats its whole analysis under all eight
consensus sequences that are equally valid under the majority vote, so no conclusion rests
on an arbitrary tie-break.

## Inputs

The raw sequencing data, the AlphaFold Server output and the intermediate count tables are
not deposited here and are available from the author on reasonable request. Each script
lists the files it expects in its docstring.

## Requirements

Python 3.10 or later, with numpy, scipy, pandas, statsmodels, openpyxl and
python-Levenshtein. The versions used were Python 3.11, NumPy 1.26, SciPy 1.11,
pandas 2.1, statsmodels 0.14 and python-Levenshtein 0.21.
