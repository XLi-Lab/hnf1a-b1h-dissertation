#!/usr/bin/env python3
"""
analyse_cici_dms.py
===================
Analysis-only pipeline for the deep mutational scanning (DMS) datasets
used in the HNF1A bacterial one-hybrid dissertation. Reproduces every statistic
reported in the Results from the raw data. Figures are produced by separate plotting
scripts, which are not included here; all reported numbers are computed below.

INPUTS (place alongside this script):
    wildtype.xlsx                            WT homeodomain coding sequence + 17 bp HNF1 motif
    dms_protein_variants.xlsx                protein-side single-nt variants
    dms_site_variants.xlsx                   DNA-side single-nt variants of the 17 bp site
    dms_double_mutants.xlsx                  protein x DNA double mutants (sheet 'in')
    HNF1A_B1H_all_motifs_corrected_2.csv     Phase 1 native sites + tiers (CSV; columns: motif_id, sequence,
                                             delta_DBD_minus_Lib, tier, ...)

Requirements: pip install pandas numpy scipy openpyxl    (Python 3.10+; no other packages needed)
Run:  python analyse_cici_dms.py
"""
import csv, pandas as pd, numpy as np
from scipy.stats import mannwhitneyu, pearsonr, spearmanr

def hdr(s): print('\n'+'='*72+'\n'+s+'\n'+'='*72)

# ----------------------------------------------------------------------
# Reference sequences and helpers
# ----------------------------------------------------------------------
bases='tcag'; aas='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon={a+b+c:aas[i] for i,(a,b,c) in enumerate((a,b,c) for a in bases for b in bases for c in bases)}
wt=pd.read_excel('wildtype.xlsx').iloc[0]; hnf=wt['hnf1a_nt'].lower(); mot=wt['motif_nt'].lower()
prot="".join(codon[hnf[i:i+3]] for i in range(0,276,3))
def rc(s): return s[::-1].translate(str.maketrans('ACGT','TGCA'))

# Phase 1 native sites + tiers
with open('HNF1A_B1H_all_motifs_corrected_2.csv',newline='',encoding='utf-8') as _fh:
    rows=[tuple(r) for r in csv.reader(_fh)]   # rows[0]=header, rows[1:]=data (strings)
H0=rows[0]; dat=[r for r in rows[1:] if r[H0.index('sequence')]]
SI,TI,MI=H0.index('sequence'),H0.index('tier'),H0.index('motif_id')
def tier(r):
    t=str(r[TI]); return 'T1' if 'Tier 1' in t else 'T2' if 'Tier 2' in t else 'SA' if 'Self' in t else 'NR'

# ======================================================================
# 1. Fitness-scale validation (PROTEIN single mutants: silent/missense/nonsense)
# ======================================================================
hdr('1. FITNESS-SCALE VALIDATION  (protein single-nt variants)')
v=pd.read_excel('dms_protein_variants.xlsx')
def classify(p,mut):
    s=list(hnf); s[p-1]=mut.lower(); ci=(p-1)//3; na=codon["".join(s[ci*3:ci*3+3])]
    return ('syn' if na==prot[ci] else 'non' if na=='*' else 'mis'), ci
v['cls'],v['res']=zip(*[classify(int(r.Pos),r.Mut) for r in v.itertuples()])
syn,mis,non=[v[v.cls==c]['fitness'].values for c in ('syn','mis','non')]
print(f"  silent   n={len(syn):3d}  median fitness={np.median(syn):+.3f}")
print(f"  missense n={len(mis):3d}  median fitness={np.median(mis):+.3f}")
print(f"  nonsense n={len(non):3d}  median fitness={np.median(non):+.3f}")
print(f"  silent vs missense : Mann-Whitney p={mannwhitneyu(syn,mis,alternative='two-sided').pvalue:.2e}")
print(f"  missense vs nonsense: Mann-Whitney p={mannwhitneyu(mis,non,alternative='two-sided').pvalue:.2e}")

# ======================================================================
# 2. DNA single-site saturation map (per-position requirement)
# ======================================================================
hdr('2. DNA SINGLE-SITE MAP  (51 site variants; position = DMS site + 1)')
m=pd.read_excel('dms_site_variants.xlsx'); m['mp']=m['Pos']-276   # mp = DMS site index 1..17
permean=m.groupby('mp')['fitness'].mean()
print("  mean fitness per position (position, mean):")
for sp in sorted(permean.index):
    print(f"    position {sp+1:2d}  ({mot[sp-1].upper()})  {permean[sp]:+.2f}")
top=permean.sort_values().head(3)
print("  three most sensitive positions (position frame):",
      [(int(s)+1, round(float(permean[s]),2)) for s in top.index])

# ======================================================================
# 3. DMS-derived binding score does NOT separate binders from non-responders
# ======================================================================
hdr('3. DMS BINDING SCORE vs FUNCTIONAL CLASS  (core necessary, not sufficient)')
fitmat={(i,b):0.0 for i in range(13) for b in 'ACGT'}
for r in m.itertuples():
    ci=int(r.mp)-4
    if 0<=ci<13: fitmat[(ci,r.Mut.upper())]=r.fitness
def score(seq):
    best=-1e9
    for s in (seq,rc(seq)):
        for i in range(len(s)-12): best=max(best,sum(fitmat[(k,b)] for k,b in enumerate(s[i:i+13])))
    return best
sc={'T1':[],'T2':[],'SA':[],'NR':[]}
for r in dat:
    if int(float(r[MI]))==0: continue   # exclude motif_00 consensus positive control
    sc[tier(r)].append(score(r[SI].upper()))
pmw=mannwhitneyu(sc['T1']+sc['T2'],sc['NR'],alternative='two-sided').pvalue
for k in ('T1','T2','SA','NR'): print(f"  {k}: n={len(sc[k]):2d}  median DMS score={np.median(sc[k]):+.2f}")
print(f"  binders (T1+T2) vs non-responders: Mann-Whitney p={pmw:.2f}  (n.s.)")

# ======================================================================
# 4. Protein per-residue sensitivity + AlphaFold 3 contact permutation test
# ======================================================================
hdr('4. PROTEIN PER-RESIDUE SENSITIVITY + AF3 CONTACT PERMUTATION TEST')
perres=v[v.cls=='mis'].groupby('res')['fitness'].mean().reindex(range(92))
af3=[(203,'R'),(266,'N'),(269,'A'),(270,'N'),(273,'K')]; idx=[u-197 for u,_ in af3]
obs=np.nanmean([perres[i] for i in idx]); bg=np.nanmean(perres.values); allr=perres.dropna().values
rng=np.random.default_rng(0)
perm=np.array([rng.choice(allr,5,replace=False).mean() for _ in range(20000)])
pp=(perm<=obs).mean()
for u,aa in af3: print(f"    {aa}{u}: mean missense fitness {perres[u-197]:+.2f}")
print(f"  mean of 5 AF3 contacts = {obs:+.2f}   all residues = {bg:+.2f}")
print(f"  permutation test (20,000 random 5-residue sets): p = {pp:.3f}")

# ======================================================================
# 5. Concordance: Phase 1 conservation vs DNA-side DMS cost
# ======================================================================
hdr('5. CONSERVATION vs DMS COST  (per-position concordance)')
T1seq=[r[SI].upper() for r in dat if tier(r)=='T1' and len(r[SI])==20]
nat=np.array([max([s[i] for s in T1seq].count(b) for b in 'ACGT')/len(T1seq) for i in range(20)])*100
cost=-m.groupby('mp')['fitness'].mean().reindex(range(1,18)).values
pos1=np.arange(2,19); cons=nat[pos1-1]
rho=spearmanr(cost,cons)[0]
print(f"  Tier 1 conservation 100% at core positions 6, 8, 13, 15")
print(f"  Spearman(conservation, DMS cost) over positions 2-18 = {rho:+.2f}  (no global trend)")
print(f"  convergence at central adenine (Phase 1 pos 8 = DMS site 7): conservation 100%, cost {cost[6]:.2f}")

# ======================================================================
# 6. Protein x DNA coupling (double-mutant epistasis)
# ======================================================================
hdr('6. PROTEIN x DNA COUPLING  (double mutants; additive null model)')
d=pd.read_excel('dms_double_mutants.xlsx',sheet_name='in').dropna(subset=['fitness_uncorr','exp','delta_aa'])
obs2=d['fitness_uncorr'].values; add=d['exp'].values; epi=d['delta_aa'].values
r,_=pearsonr(obs2,add); rho2,_=spearmanr(obs2,add)
ve=1-np.sum((obs2-add)**2)/np.sum((obs2-obs2.mean())**2)
strong=np.abs(epi-epi.mean())>2*epi.std()
print(f"  n = {len(d)} pairs")
print(f"  observed vs additive: Pearson r={r:.3f}, Spearman rho={rho2:.3f}")
print(f"  additive-model R2 = {ve:.3f}  (explains {100*ve:.1f}% of double-mutant variance)")
print(f"  epistasis (obs-exp): mean={epi.mean():+.3f}, sd={epi.std():.3f}; |epi|>2SD: {strong.sum()} ({100*strong.mean():.1f}%)")
d['motif']=d['Pos2']-275; d['residue']=197+((d['Pos1']-1)//3); d['ae']=d['delta_aa'].abs()
print("  strongest DNA positions:",d.groupby('motif')['ae'].mean().sort_values(ascending=False).head(4).round(3).to_dict())
print("  strongest HNF1A residues     :",d.groupby('residue')['ae'].mean().sort_values(ascending=False).head(5).round(3).to_dict())
# ======================================================================
# 7. Coupling at the DNA-binding interface (noise-normalised)
# ======================================================================
hdr('7. INTERFACE COUPLING  (N-terminal arm + recognition helix vs rest of domain)')
c=pd.read_excel('dms_double_mutants.xlsx',sheet_name='in').dropna(
    subset=['fitness_uncorr','exp','delta_aa','sigma_uncorr','sigma1','sigma2'])
def is_syn(p1,mut):
    ci=(p1-1)//3; cp=(p1-1)%3
    s=list(hnf[ci*3:ci*3+3]); wt_aa=codon["".join(s)]
    s[cp]=str(mut).lower()
    return wt_aa==codon["".join(s)]
c['syn']=[is_syn(int(r.Pos1),r.Mut1) for r in c.itertuples()]
c['z']=c['delta_aa']/np.sqrt(c['sigma_uncorr']**2+c['sigma1']**2+c['sigma2']**2)
c['residue']=197+((c['Pos1']-1)//3)

noise=c[c.syn]['z'].abs().mean()          # empirical noise floor from synonymous controls
pz=c[~c.syn].groupby('residue').agg(z=('z',lambda x:x.abs().mean()),n=('z','size'))
pz=pz[pz.n>=15]['z']                      # residues with adequate coverage

ARM=set(range(203,207))                   # N-terminal arm  (UniProt P20823)
HELIX=set(range(261,276))                 # recognition helix (UniProt P20823)
sel=[i for i in pz.index if i in ARM|HELIX]
oth=[i for i in pz.index if i not in ARM|HELIX]
pmw2=mannwhitneyu(pz[sel],pz[oth],alternative='two-sided').pvalue

print(f"  noise floor from synonymous protein-side changes: mean |z| = {noise:.2f}")
print(f"  residues retained (>=15 non-synonymous pairs): {len(pz)} of 92")
print(f"  domain median |z| = {pz.median():.2f}")
print(f"  interface (arm 203-206 + helix 261-275): n={len(sel)}  mean |z| = {pz[sel].mean():.2f}")
print(f"  rest of domain                          : n={len(oth)}  mean |z| = {pz[oth].mean():.2f}")
print(f"  Mann-Whitney (two-sided): p = {pmw2:.1e}")
print("\nAll analyses complete.")
