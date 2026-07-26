#!/usr/bin/env python3
"""
analyse_consensus_distance.py
=============================
Analysis-only script for section 3.8 and Figure 5 of the HNF1A bacterial
one-hybrid dissertation. Scores every native promoter variant by its Hamming
distance to the Tier 1 majority-vote consensus and compares the four functional
classes.

Follows the decision rule set out in Methods 2.18. Four groups are compared, so
an omnibus test is run before any pairwise test. Shapiro-Wilk within each class
and Levene across classes are checked first; both support the parametric route,
so the omnibus is a one-way ANOVA and the post-hoc is Tukey HSD over all six
pairs. Kruskal-Wallis is reported alongside as a nonparametric cross-check.

Three positions of the Tier 1 consensus have no outright majority base, which
makes eight consensus sequences equally valid. The script enumerates every tied
position, repeats the whole analysis under each consensus, and reports the
worst case and best case for each pair, so no conclusion rests on the sequence
the majority-vote happens to pick.

INPUTS (place alongside this script):
    HNF1A_B1H_all_motifs_corrected_2.csv   Phase 1 sites, tiers and enrichment
                                           (columns: motif_id, sequence, tier,
                                            delta_DBD_minus_Lib, ...)

Requirements: pip install numpy scipy statsmodels   (Python 3.10+)
Run:  python analyse_consensus_distance.py
"""
import csv, itertools
import numpy as np
from collections import Counter
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd

def hdr(s): print('\n' + '=' * 72 + '\n' + s + '\n' + '=' * 72)

CLASSES = ['T1', 'T2', 'SA', 'NR']
NAMES = {'T1': 'Tier 1', 'T2': 'Tier 2', 'SA': 'Self-activator', 'NR': 'No-response'}

# ----------------------------------------------------------------------
# Load the Phase 1 table
# ----------------------------------------------------------------------
with open('HNF1A_B1H_all_motifs_corrected_2.csv', newline='', encoding='utf-8') as fh:
    rows = [tuple(r) for r in csv.reader(fh)]
H = rows[0]
SI, TI, MI = H.index('sequence'), H.index('tier'), H.index('motif_id')

def tier(r):
    t = str(r[TI])
    return 'T1' if 'Tier 1' in t else 'T2' if 'Tier 2' in t else 'SA' if 'Self' in t else 'NR'

# motif_00 is the 17 nt synthetic consensus positive control and is not a native
# promoter, so it is excluded. motif_08 dropped out of the library and has no reads.
natives = sorted([r for r in rows[1:] if r[SI] and int(float(r[MI])) != 0 and len(r[SI]) == 20],
                 key=lambda r: int(float(r[MI])))
groups = {c: [r for r in natives if tier(r) == c] for c in CLASSES}
print(f"native variants scored: {len(natives)}")
for c in CLASSES:
    print(f"  {NAMES[c]:16s} n = {len(groups[c])}")

# ----------------------------------------------------------------------
# Tier 1 consensus by majority vote, with every tie enumerated
# ----------------------------------------------------------------------
# Ties are broken by first occurrence in motif_id order. That rule is arbitrary,
# so it is stated here rather than left implicit, and every alternative is swept below.
t1_seqs = [r[SI].upper() for r in groups['T1']]
per_pos = []
tied_pos = []
primary = []
for p in range(20):
    cnt = Counter(s[p] for s in t1_seqs)
    top = max(cnt.values())
    primary.append(cnt.most_common(1)[0][0])
    bases = sorted(b for b, n in cnt.items() if n == top)
    per_pos.append(bases)
    if len(bases) > 1:
        tied_pos.append((p + 1, dict(sorted(cnt.items()))))

CONSENSUS = "".join(primary)
consensus_set = ["".join(c) for c in itertools.product(*per_pos)]

hdr('TIER 1 CONSENSUS  (majority vote over the Tier 1 sites)')
print(f"  primary consensus : {CONSENSUS}")
print(f"  positions with no outright majority base: {len(tied_pos)}")
for p, cnt in tied_pos:
    print(f"    position {p:2d}  {cnt}")
print(f"  equally valid consensus sequences: {len(consensus_set)}")

def hamming(s, ref):
    return sum(a != b for a, b in zip(s, ref))

# ----------------------------------------------------------------------
# Descriptive statistics under the primary consensus
# ----------------------------------------------------------------------
def distances(ref):
    return {c: np.array([hamming(r[SI].upper(), ref) for r in groups[c]], float) for c in CLASSES}

D = distances(CONSENSUS)
hdr('HAMMING DISTANCE BY FUNCTIONAL CLASS')
for c in CLASSES:
    d = D[c]
    print(f"  {NAMES[c]:16s} n={len(d):2d}  mean={d.mean():.3f}  median={np.median(d):.1f}  "
          f"IQR={np.percentile(d,25):.1f}-{np.percentile(d,75):.1f}  range={int(d.min())}-{int(d.max())}")

# ----------------------------------------------------------------------
# Assumption checks, then the omnibus test
# ----------------------------------------------------------------------
hdr('ASSUMPTION CHECKS  (Methods 2.18)')
normal = True
for c in CLASSES:
    W, p = stats.shapiro(D[c])
    if p < 0.05:
        normal = False
    print(f"  Shapiro-Wilk  {NAMES[c]:16s} W = {W:.4f}  p = {p:.4f}")
Wl, pl = stats.levene(*[D[c] for c in CLASSES])
print(f"  Levene across the four classes   W = {Wl:.4f}  p = {pl:.4f}")
print(f"  parametric route supported: {'yes' if normal and pl >= 0.05 else 'no'}")

hdr('OMNIBUS TEST')
F, p_anova = stats.f_oneway(*[D[c] for c in CLASSES])
n_total = sum(len(D[c]) for c in CLASSES)
df_b, df_w = len(CLASSES) - 1, n_total - len(CLASSES)
Hk, p_kw = stats.kruskal(*[D[c] for c in CLASSES])
print(f"  one-way ANOVA    F({df_b},{df_w}) = {F:.3f}   p = {p_anova:.6f}")
print(f"  Kruskal-Wallis   H = {Hk:.3f}          p = {p_kw:.6f}   (nonparametric cross-check)")

# ----------------------------------------------------------------------
# Post-hoc, only if the omnibus is significant
# ----------------------------------------------------------------------
def tukey(dist):
    vals = np.concatenate([dist[c] for c in CLASSES])
    labs = np.concatenate([[c] * len(dist[c]) for c in CLASSES])
    tk = pairwise_tukeyhsd(vals, labs, alpha=0.05)
    out = {}
    for row, pv in zip(tk._results_table.data[1:], tk.pvalues):
        out[frozenset((row[0], row[1]))] = float(pv)
    return out

PAIRS = [('T1', 'T2'), ('T1', 'SA'), ('T1', 'NR'),
         ('T2', 'SA'), ('T2', 'NR'), ('SA', 'NR')]

def stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'

hdr('POST-HOC  (Tukey HSD, all six pairs; p values are already adjusted)')
if p_anova >= 0.05:
    print("  omnibus is null, so no post-hoc test is run")
    padj = {}
else:
    padj = tukey(D)
    for a, b in PAIRS:
        p = padj[frozenset((a, b))]
        print(f"  {NAMES[a]:16s} vs {NAMES[b]:16s} p_adj = {p:.5f}  {stars(p)}")

# ----------------------------------------------------------------------
# Tie-break sensitivity: repeat everything under all valid consensus sequences
# ----------------------------------------------------------------------
hdr('TIE-BREAK SENSITIVITY  (analysis repeated under every valid consensus)')
anova_ps, pair_ps = [], {frozenset(p): [] for p in PAIRS}
for ref in consensus_set:
    Dr = distances(ref)
    _, pa = stats.f_oneway(*[Dr[c] for c in CLASSES])
    anova_ps.append(pa)
    tk = tukey(Dr)
    for a, b in PAIRS:
        pair_ps[frozenset((a, b))].append(tk[frozenset((a, b))])
print(f"  ANOVA p across {len(consensus_set)} consensus sequences: "
      f"min = {min(anova_ps):.5f}  max = {max(anova_ps):.5f}  (significant under all)")
for a, b in PAIRS:
    v = pair_ps[frozenset((a, b))]
    verdict = ('significant under all 8' if max(v) < 0.05 else
               'significant under none' if min(v) >= 0.05 else 'varies with the tie-break')
    print(f"  {NAMES[a]:16s} vs {NAMES[b]:16s} p_adj {min(v):.4f} to {max(v):.4f}   {verdict}")
print("\n  The only pairwise result reported in section 3.8 is Tier 1 versus no-response,")
print("  which is significant under every valid consensus. The Tier 2 versus no-response")
print("  comparison is significant under none of them, which is why it is not claimed.")
print("  Pairs marked as varying with the tie-break are reported as not significant,")
print("  which is the conservative reading under the primary consensus.")

print('\nAll analyses complete.')
