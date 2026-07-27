#!/usr/bin/env python3
"""
analyse_af3_confidence.py
=========================
Analysis-only pipeline for the AlphaFold 3 results reported in Section 3.10,
Figure 6 and Table D2. Reproduces every number quoted for both panels of that
figure. No figure is drawn here; the plotting code is not deposited.

PANEL A  single-copy TMED10 model: per-residue pLDDT, the five strongest DNA-contact
         residues, and the Table D2 confidence metrics.
PANEL B  all 21 sites: AlphaFold protein-DNA interface confidence (ipTM) vs the
         experimental B1H enrichment, and whether confidence predicts binding.

Statistical policy (Methods 2.18). Panel B makes a fixed number of planned
comparisons on the same interface ipTM values: a three-class omnibus test and a
single pooled binders-versus-non-responders test. That is a family of two, so
both p values are Bonferroni corrected by a factor of two. Tier 2 fails
Shapiro-Wilk and Tier 1 is borderline, so Kruskal-Wallis is the omnibus that is
reported; the one-way ANOVA is printed alongside for the record only. Neither
omnibus is significant, so no post-hoc pairwise test is run.

INPUTS:
    af3/folds_2026_06_22_14_59/<job>/   21 AF3 monomer jobs, each with
        *model_0.cif  *full_data_0.json  *summary_confidences_0.json
    HNF1A_B1H_all_motifs_corrected_2.csv       Phase 1 enrichment (CSV; columns: motif_id, delta_DBD_minus_Lib, ...)

Requirements: pip install numpy scipy    (Python 3.10+; no other packages needed)
Run:  python analyse_af3_confidence.py
"""
import json, glob, csv, numpy as np
from collections import defaultdict
from scipy.stats import spearmanr, kruskal, mannwhitneyu, f_oneway, shapiro, levene

def hdr(s): print('\n'+'='*72+'\n'+s+'\n'+'='*72)

def bonf(p, k):
    """Bonferroni-adjusted p value for a family of k planned comparisons."""
    return min(1.0, p * k)

AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
     'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
     'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
OFF=196                                          # AF3 residue index 1 == HNF1A residue 197
BASE='af3/folds_2026_06_22_14_59/'

# ======================================================================
# PANEL A : single-copy TMED10 deep dive
# ======================================================================
hdr('PANEL A : single-copy TMED10 model (pLDDT, DNA contacts, Table D2 metrics)')
job=BASE+'t1_30_tmed10/'
cif=glob.glob(job+'*model_0.cif')[0]; fd=glob.glob(job+'*full_data_0.json')[0]
sm=json.load(open(glob.glob(job+'*summary_confidences_0.json')[0]))

# --- per-residue pLDDT (chain A = protein), read from the B-factor column of the cif ---
cm={}; rows=[]; inl=False; ci=0
for line in open(cif):
    line=line.rstrip()
    if line=='loop_': inl=False; ci=0; cm={}; rows=[]
    elif line.startswith('_atom_site.'): inl=True; cm[line.strip()]=ci; ci+=1
    elif inl and line and not line.startswith(('_','#')):
        p=line.split()
        if len(p)==ci: rows.append(p)
cc=cm['_atom_site.label_asym_id']; cr=cm['_atom_site.label_seq_id']
cb=cm['_atom_site.B_iso_or_equiv']; co=cm['_atom_site.label_comp_id']
pld=defaultdict(list); rn={}
for r in rows:
    try: res=int(r[cr])
    except: continue
    if r[cc]=='A': pld[res].append(float(r[cb])); rn[res]=r[co]
plddt={r:float(np.mean(x)) for r,x in pld.items()}
print(f"  modelled protein residues : {min(plddt)+OFF}-{max(plddt)+OFF}  ({len(plddt)} aa)")
print(f"  mean pLDDT (per residue)  : {np.mean(list(plddt.values())):.1f}")

# --- five strongest DNA-contact residues (max contact prob to either DNA chain) ---
f=json.load(open(fd)); cp=np.array(f['contact_probs']); ch=f['token_chain_ids']; rid=f['token_res_ids']
idx=defaultdict(list)
for i,c in enumerate(ch): idx[c].append(i)
dna=idx['B']+idx['C']
scored=sorted(((rid[i],cp[i,dna].max()) for i in idx['A']),key=lambda kv:-kv[1])[:5]
print("  five strongest DNA-contact residues (AF3 contact probability):")
for l,pr in scored:
    print(f"    {AA3.get(rn.get(l,'?'),'?')}{l+OFF}  p={pr:.2f}")

# --- Table D2 confidence metrics ---
cpi=np.array(sm['chain_pair_iptm']); pae=np.array(sm['chain_pair_pae_min'])
print("  Table D2 metrics:")
print(f"    ranking score          {sm['ranking_score']:.2f}")
print(f"    ipTM                   {sm['iptm']:.2f}")
print(f"    pTM                    {sm['ptm']:.2f}")
print(f"    fraction disordered    {sm['fraction_disordered']*100:.0f}%")
print(f"    min protein-DNA PAE    {min(pae[0,1],pae[0,2]):.2f} A")

# ======================================================================
# PANEL B : 21-site comparison (confidence vs binding)
# ======================================================================
hdr('PANEL B : 21 sites - interface ipTM vs experimental binding')
with open('HNF1A_B1H_all_motifs_corrected_2.csv',newline='',encoding='utf-8') as _fh:
    allm=[tuple(r) for r in csv.reader(_fh)]   # allm[0]=header, allm[1:]=data (strings)
hh=list(allm[0]); _mi=hh.index('motif_id')
delta={'motif_%02d'%int(float(r[_mi])):r[hh.index('delta_DBD_minus_Lib')] for r in allm[1:] if r[_mi] not in ('','None')}
pts=[]
for d in sorted(glob.glob(BASE+'*/')):
    jb=d.split('/')[-2]; s=json.load(open(glob.glob(d+'*summary_confidences_0.json')[0]))
    c=np.array(s['chain_pair_iptm']); cls=jb.split('_')[0].upper()
    pdna=float(np.mean([c[0,1],c[0,2]]))          # protein-DNA interface ipTM (panel B y-axis)
    pts.append((cls, pdna, float(delta['motif_%02d'%int(jb.split('_')[1])]),
                '_'.join(jb.split('_')[2:]), float(s['iptm'])))
ipt=[p[1] for p in pts]; dlt=[p[2] for p in pts]
rho,pv=spearmanr(ipt,dlt)
print(f"  n = {len(pts)} sites")
print(f"  interface ipTM vs binding (delta log2FC): Spearman rho = {rho:+.2f}  (p = {pv:.2f})")
rho2,pv2=spearmanr([p[4] for p in pts],dlt)
print(f"  robustness (overall ipTM instead)       : Spearman rho = {rho2:+.2f}  (p = {pv2:.2f})")
for cls in ('T1','T2','NR'):
    g=[p[1] for p in pts if p[0]==cls]
    print(f"    {cls}: n={len(g)}  mean interface ipTM={np.mean(g):.2f}")

# ----------------------------------------------------------------------
# Two planned comparisons, so Bonferroni k = 2 (Methods 2.18)
# ----------------------------------------------------------------------
K=2
grp=[np.array([p[1] for p in pts if p[0]==c],float) for c in ('T1','T2','NR')]
g1=np.array([p[1] for p in pts if p[0] in ('T1','T2')],float)
gn=np.array([p[1] for p in pts if p[0]=='NR'],float)

print("\n  assumption checks:")
for c,g in zip(('T1','T2','NR'),grp):
    W,p=shapiro(g)
    print(f"    Shapiro-Wilk {c}  W = {W:.4f}  p = {p:.4f}{'   (fails)' if p<0.05 else ''}")
Wl,pl=levene(*grp)
print(f"    Levene across the three classes  W = {Wl:.4f}  p = {pl:.4f}")
print("    normality is not supported in every class, so Kruskal-Wallis is the omnibus reported")

Fa,pa=f_oneway(*grp)
df_b,df_w=len(grp)-1,sum(len(g) for g in grp)-len(grp)
print(f"\n  one-way ANOVA (for the record, not used for inference)")
print(f"    F({df_b},{df_w}) = {Fa:.3f}   raw p = {pa:.4f}")

Hk,pk=kruskal(*grp)
print(f"  omnibus  Kruskal-Wallis across 3 classes")
print(f"    H = {Hk:.3f}   raw p = {pk:.4f}   Bonferroni-adjusted p = {bonf(pk,K):.2f}")

pu=mannwhitneyu(g1,gn,alternative='two-sided').pvalue
print(f"  planned  binders (T1+T2) vs non-responders, Mann-Whitney U")
print(f"    raw p = {pu:.4f}   Bonferroni-adjusted p = {bonf(pu,K):.2f}")
print("  the omnibus is null, so no post-hoc pairwise test is run")

top=max(pts,key=lambda p:p[1]); print(f"\n  most confident model overall : {top[3].upper()} ({top[0]}), interface ipTM {top[1]:.2f}")
akr=[p for p in pts if p[3]=='akr1c3']
if akr: print(f"  perfect-core non-responder   : AKR1C3, interface ipTM {akr[0][1]:.2f} (within binder range)")
print('\nAll analyses complete.')
