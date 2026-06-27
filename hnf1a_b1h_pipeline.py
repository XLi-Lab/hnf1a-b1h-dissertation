#!/usr/bin/env python3
"""
HNF1A B1H Phase 1: full analysis pipeline

Takes raw Nanopore amplicon FASTQs through to the tissue-annotated
3-panel figure. Demultiplexing follows the original analysis, so the
counts reproduce.

Inputs (place in ./input/):
    Lib-only_pool_Tube1_raw_fastq.gz     (FCL order #27683)
    DBD-Lib_pool_Tube2_raw_fastq.gz      (FCL order #27683)
    hnf1a_63_unique_motif.txt            (63 designed library variants)
    blat_best_hits.json                  (BLAT results vs hg38)
    BLAT_gene_overlaps.tsv               (UCSC Table Browser output)

Outputs (in ./output/):
    HNF1A_B1H_3panel_VERIFIED.png
    HNF1A_B1H_VERIFIED_ordered_with_hits.tsv
    HNF1A_B1H_VERIFIED_annotations.tsv

Dependencies:
    python >= 3.9
    numpy, matplotlib, pandas, python-Levenshtein
"""

import os
import sys
import gzip
import math
import json
import csv
import re
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

try:
    import Levenshtein as _lev
    def levenshtein(a, b):
        return _lev.distance(a, b)
except ImportError:
    print("WARNING: python-Levenshtein not installed. Using pure-Python fallback (slower).")
    print("         Install with: pip install python-Levenshtein")
    def levenshtein(a, b):
        if abs(len(a) - len(b)) > 10:
            return abs(len(a) - len(b))
        prev = list(range(len(b) + 1))
        for i, ca in enumerate(a, 1):
            cur = [i]
            for j, cb in enumerate(b, 1):
                cost = 0 if ca == cb else 1
                cur.append(min(cur[-1]+1, prev[j]+1, prev[j-1]+cost))
            prev = cur
        return prev[-1]


INPUT_DIR = Path("./input")
OUTPUT_DIR = Path("./output")
OUTPUT_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Sample-to-barcode mapping (lab notebook section 28.3)
# ---------------------------------------------------------------------------
# Pool 1 (Tube1) = Lib-only  |  Pool 2 (Tube2) = DBD-Lib
# Inputs are unbarcoded. They are detected by the forward primer prefix
# "CCGCTCAT" (the 5' end of the amplicon that has no 8 nt sample barcode).

SAMPLE_BARCODES = {
    "Lib_input": None,       "Lib1": "GACACACT", "Lib2": "ACGACTTG", "Lib3": "ATCGATCG",
    "DBD_input": None,       "TF1":  "GCCGAATT", "TF2":  "CGGCAATT", "TF3":  "GAACGTTC",
}
SAMPLE_POOL = {
    "Lib_input": "Tube1", "Lib1": "Tube1", "Lib2": "Tube1", "Lib3": "Tube1",
    "DBD_input": "Tube2", "TF1":  "Tube2", "TF2":  "Tube2", "TF3":  "Tube2",
}
FORWARD_PRIMER_5_PREFIX = "CCGCTCAT"

FLANK_5 = "GCGGCCGC"  # NotI
FLANK_3 = "GAATTC"    # EcoRI

T1_DELTA = 2.5
T2_DELTA = 1.5
PSEUDO = 1.0
POS_CONTROL_SEQ = "AGTTAATTATTAACCAA"


# ---------------------------------------------------------------------------
# Odom 2004 Supplementary Tables S2 + S3 (embedded)
# ---------------------------------------------------------------------------
# Each tuple is (gene_symbol_as_printed_in_Odom_2004, RefSeq_accession).
# Used for cross-referencing our B1H hits against known HNF1a ChIP targets
# in primary hepatocytes (S2) and pancreatic islets (S3).

ODOM_S2_HEPATOCYTE = [
    ("AADAC","NM_001086"),("ABCC2","NM_000392"),("ACF","NM_014576"),("ADH1A","NM_000667"),
    ("ADH1B","NM_000668"),("ADH6","NM_000672"),("AGT","NM_000029"),("AHSG","NM_001622"),
    ("AK2","NM_001625"),("AKR1C2","NM_001354"),("AKR1C3","NM_003739"),("AKR1C4","NM_001818"),
    ("ALB","NM_000477"),("ALDH3A2","NM_000382"),("ALS2","NM_020919"),("AMBP","NM_001633"),
    ("ANGPTL3","NM_014495"),("ANPEP","NM_001150"),("AP3M1","NM_012095"),("APCS","NM_001639"),
    ("APG3","NM_022488"),("APOA2","NM_001643"),("APOH","NM_000042"),("AQP3","NM_004925"),
    ("AQP9","NM_020980"),("ARHGAP11A","NM_014783"),("ASGR1","NM_001671"),("ASGR2","NM_001181"),
    ("ATF2","NM_001880"),("AUTL1","NM_032852"),("BAT3","NM_004639"),("BIKE","NM_017593"),
    ("BTN2A1","NM_078476"),("C1S","NM_001734"),("C2","NM_000063"),("C4BPA","NM_000715"),
    ("C8B","NM_000066"),("CCNE1","NM_001238"),("CDCA1","NM_031423"),("CISH","NM_013324"),
    ("CLYBL","NM_138280"),("CNTNAP2","NM_014141"),("CPB2","NM_016413"),("CREBL2","NM_001310"),
    ("CRP","NM_000567"),("CTSZ","NM_001336"),("CYB5","NM_001914"),("CYB5-M","NM_030579"),
    ("CYP2E","NM_000773"),("CYP3A43","NM_022820"),("DAF","NM_000574"),("DC13","NM_020188"),
    ("DKFZP564O0463","NM_014156"),("DKFZP586A0522","NM_014033"),("DKFZP586M0122","NM_015425"),
    ("DLEU1","NM_005887"),("DUSP6","NM_022652"),("EIF4EBP2","NM_004096"),("ELF3","NM_004433"),
    ("ENPEP","NM_001977"),("F11","NM_019559"),("FE65L2","NM_006051"),("FH","NM_000143"),
    ("FKSG87","NM_032029"),("FLJ10242","NM_018036"),("FLJ10276","NM_018045"),("FLJ10525","NM_018126"),
    ("FLJ10583","NM_018148"),("FLJ10650","NM_018168"),("FLJ10774","NM_024662"),("FLJ11000","NM_018295"),
    ("FLJ11838","NM_024664"),("FLJ12788","NM_022492"),("FLJ13448","NM_025147"),("FLJ13611","NM_024941"),
    ("FLJ14356","NM_030824"),("FLJ20080","NM_017657"),("FLJ20718","NM_017939"),("FLJ21272","NM_025032"),
    ("FLJ21934","NM_024743"),("FLJ22551","NM_024708"),("FLJ23259","NM_024727"),("FNTB","NM_002028"),
    ("G0S2","NM_015714"),("G3A","NM_019101"),("G6PT1","NM_001467"),("GARS","NM_002047"),
    ("GBE1","NM_000158"),("GCKR","NM_001486"),("GDI2","NM_001494"),("GIOT-2","NM_016264"),
    ("GJB1","NM_000166"),("GOT1","NM_002079"),("GPR39","NM_001508"),("GPX2","NM_002083"),
    ("GRHPR","NM_012203"),("GTF2B","NM_001514"),("GTF2E1","NM_005513"),("GTPBG3","NM_032620"),
    ("HABP2","NM_004132"),("HAL","NM_002108"),("HAO1","NM_017545"),("HCAP-G","NM_022346"),
    ("HGD","NM_000187"),("HGFAC","NM_001528"),("HNF4A","NM_000457"),("HNF4a7","AF509467"),
    ("HNMT","NM_006895"),("HPCL2","NM_012260"),("HPX","NM_000613"),("HSD11B1","NM_005525"),
    ("HSD17B2","NM_002153"),("HSPC111","NM_016391"),("HSPC129","NM_016396"),("IFNAR1","NM_000629"),
    ("IGF1R","NM_000875"),("IGFBP1","NM_000596"),("INADL","NM_005799"),("ITIH3","NM_002217"),
    ("ITIH4","NM_002218"),("ITM2B","NM_021999"),("KIAA0022","NM_014880"),("KIAA0669","NM_014779"),
    ("KIAA0844","NM_014951"),("KIAA0872","NM_014940"),("KIAA1041","NM_014947"),("KNG","NM_000893"),
    ("LBP","NM_004139"),("LOC51060","NM_015913"),("LOC51096","NM_016001"),("LOC51326","NM_016632"),
    ("LOC54518","NM_019043"),("LOC56902","NM_020143"),("LOC58486","NM_021211"),("LY6E","NM_002346"),
    ("M17S2","NM_031858"),("M96","NM_007358"),("MAGEA9","NM_005365"),("MGC10500","NM_031477"),
    ("MGC11034","NM_031453"),("MGC11266","NM_024322"),("MGC13010","NM_032687"),("MGC15435","NM_032367"),
    ("MGC955","NM_024097"),("MIA2","NM_054024"),("MRPL15","NM_014175"),("MRPS18B","NM_014046"),
    ("MSH6","NM_000179"),("MT1H","NM_005951"),("MT1L","NM_002450"),("MT1X","NM_005952"),
    ("MTHFD1","NM_005956"),("MTP","NM_000253"),("NAPA","NM_003827"),("NET-2","NM_012338"),
    ("NFKBIB","NM_002503"),("NPC1L1","NM_013389"),("NR0B2","NM_021969"),("NR1D1","NM_021724"),
    ("NR5A2","NM_003822"),("NRD1","NM_002525"),("PAFAH2","NM_000437"),("PAX8","NM_013952"),
    ("PCK1","NM_002591"),("PHF2","NM_005392"),("PIST","NM_020399"),("PLCB1","NM_015192"),
    ("PLG","NM_000301"),("PLGL","NM_002665"),("PS-PLA1","NM_015900"),("PZP","NM_002864"),
    ("RAB33B","NM_031296"),("RAMP","NM_016448"),("RARB","NM_016152"),("RBP5","NM_031491"),
    ("RNGTT","NM_003800"),("RPL37AP1","NG_000988"),("SAC","NM_018417"),("SCYE1","NM_004757"),
    ("SEL1L","NM_005065"),("SERPINA1","NM_000295"),("SERPINA10","NM_016186"),("SERPINA6","NM_001756"),
    ("SERPINC1","NM_000488"),("SERPINE1","NM_000602"),("SERPING1","NM_000062"),("SGK2","NM_016276"),
    ("SLC17A2","NM_005835"),("SLC22A11","NM_018484"),("SLPI","NM_003064"),("SNX17","NM_014748"),
    ("SRI","NM_003130"),("SSA2","NM_004600"),("SSTR1","NM_001049"),("SSTR4","NM_001052"),
    ("STRAIT11499","NM_021242"),("SUPV3L1","NM_003171"),("SYN3","NM_133632"),("TARS","NM_003191"),
    ("TBPL1","NM_004865"),("TEF","NM_003216"),("TFRC","NM_003234"),("TIEG2","NM_003597"),
    ("TM4SF4","NM_004617"),("TMEM1","NM_003274"),("TNFRSF6","NM_000043"),("UGT1A1","NM_000463"),
    ("UGT2B11","NM_001073"),("UGT2B15","NM_001076"),("UQCRC2","NM_003366"),("VNN3","NM_018399"),
    ("VTN","NM_000638"),("WBP4","NM_007187"),("WDF2","NM_052950"),("WDR12","NM_018256"),
    ("XDH","NM_000379"),("XPC","NM_004628"),("ZK1","NM_005815"),("ZNF288","NM_015642"),
    ("ZNF361","NM_018555"),
]

ODOM_S3_ISLET = [
    ("AADAC","NM_001086"),("ABCC9","NM_020297"),("ADH4","NM_000670"),("APOH","NM_000042"),
    ("ARHGAP11A","NM_014783"),("B29","NM_031939"),("BCL6","NM_001706"),("BIKE","NM_017593"),
    ("C4BPA","NM_000715"),("C6orf11","NM_005452"),("CDC45L","NM_003504"),("COL3A1","NM_000090"),
    ("COQ7","NM_016138"),("CPXCR1","NM_033048"),("CRH","NM_000756"),("CTSZ","NM_001336"),
    ("CYB5-M","NM_030579"),("DKFZP564J157","NM_018457"),("DLEU1","NM_005887"),("DOCK1","NM_001380"),
    ("DSC1","NM_024421"),("EIF3S6","NM_001568"),("ELF3","NM_004433"),("FBXO8","NM_012180"),
    ("FE65L2","NM_006051"),("FIL1","NM_014440"),("FLJ10242","NM_018036"),("FLJ10252","NM_018040"),
    ("FLJ10474","NM_018104"),("FLJ10650","NM_018168"),("FLJ11301","NM_018385"),("FLJ13273","NM_024751"),
    ("FLJ13385","NM_024853"),("FLJ13448","NM_025147"),("FLJ14855","NM_033210"),("FLJ20156","NM_017691"),
    ("FLJ20225","NM_019062"),("FLJ20234","NM_017720"),("FLJ20298","NM_017752"),("FLJ20643","NM_017916"),
    ("FLJ20731","NM_017946"),("FLJ21272","NM_025032"),("FLJ22559","NM_024928"),("FNTB","NM_002028"),
    ("GCNT3","NM_004751"),("GIOT-2","NM_016264"),("GLA","NM_000169"),("GNB2L1","NM_006098"),
    ("GPR74","NM_004885"),("H4F2","NM_003548"),("HAVCR-1","NM_012206"),("HHLA2","NM_007072"),
    ("HNF4a7","AF509467"),("IFNA10","NM_002171"),("INSR","NM_000208"),("KIAA0101","NM_014736"),
    ("KIAA0399","NM_015113"),("KIAA0844","NM_014951"),("KIF13A","NM_022113"),("KIR-023GB","NM_015868"),
    ("KIR2DS2","NM_012312"),("KIR3DL1","NM_013289"),("KRTAP1.1","NM_030967"),("KRTHA3A","NM_004138"),
    ("LIPA","NM_000235"),("LOC113201","NM_138423"),("LOC113220","NM_138424"),("LOC51092","NM_015996"),
    ("LOC56906","NM_020147"),("MCCC1","NM_020166"),("MGC10500","NM_031477"),("MGC15677","NM_032878"),
    ("MIA2","NM_054024"),("MRPL15","NM_014175"),("NPY2R","NM_000910"),("NR0B2","NM_021969"),
    ("NR2C2","NM_003298"),("NR5A2","NM_003822"),("PAFAH2","NM_000437"),("PAX8","NM_013952"),
    ("PEX13","NM_002618"),("PGCP","NM_016134"),("PRO2032","NM_018615"),("PSMA5","NM_002790"),
    ("PS-PLA1","NM_015900"),("RAB33B","NM_031296"),("RAB6KIFL","NM_005733"),("SDCCAG10","NM_005869"),
    ("SEL1L","NM_005065"),("SGK2","NM_016276"),("SLC26A7","NM_052832"),("SPO11","NM_012444"),
    ("SRI","NM_003130"),("SSTR1","NM_001049"),("TACR3","NM_001059"),("TM4SF4","NM_004617"),
    ("TMOD2","NM_014548"),("TMP21","NM_006827"),("UQCRC2","NM_003366"),("UROD","NM_000374"),
    ("VNN3","NM_018399"),("WBP4","NM_007187"),("ZNF155","NM_003445"),("ZNF300","NM_052860"),
    ("pcnp","NM_020357"),("Nod1(-)6kb","NM_006092"),
]


# ---------------------------------------------------------------------------
# STAGE 1: Read FASTQ, demultiplex, map reads to library variants
# ---------------------------------------------------------------------------
# Demultiplexing logic:
#  - orient each read so NotI comes before EcoRI (try both strands)
#  - read the 8 nt barcode at offset 0 or 1 (allows a 1 bp frame shift)
#  - allow Hamming distance <= 1 (tolerates one barcode sequencing error)
#  - a Tube1 read is only matched against Tube1 barcodes, and vice versa
#  - unbarcoded input reads are detected by the forward primer prefix

def load_library(path):
    """Load the 63 designed library motifs from the CSV file, plus motif_0 as positive control."""
    lib = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(",")
            if len(parts) >= 2:
                idx = int(parts[0])
                lib[idx] = parts[1].strip()
    lib[0] = POS_CONTROL_SEQ
    return lib

def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTN","TGCAN"))[::-1]

def hamming_distance(a, b):
    return sum(x != y for x, y in zip(a, b))

def get_canonical_orientation(read):
    """Return (canonical_read, variable_region) where NotI precedes EcoRI.
    Tries both read orientations. Variable region must be 14-23 bp."""
    for candidate in (read, reverse_complement(read)):
        i = candidate.find(FLANK_5)
        if i < 0:
            continue
        j = candidate.find(FLANK_3, i + len(FLANK_5))
        if j < 0:
            continue
        var = candidate[i + len(FLANK_5) : j]
        if 14 <= len(var) <= 23:
            return candidate, var
    return None, None

def identify_sample(canonical_read, pool, max_mismatches=1):
    """Read the 8 nt barcode at offsets 0 or 1 from the canonical read.
    Only considers barcodes belonging to the given pool."""
    candidates = [s for s in SAMPLE_BARCODES if SAMPLE_POOL[s] == pool]
    for offset in (0, 1):
        prefix = canonical_read[offset:offset + 8]
        if len(prefix) < 8:
            continue
        for sample in candidates:
            barcode = SAMPLE_BARCODES[sample] or FORWARD_PRIMER_5_PREFIX
            if hamming_distance(prefix, barcode) <= max_mismatches:
                return sample
    return None

def assign_to_motif(var_region, lib_seqs, max_edit=2):
    best_id, best_dist = None, max_edit + 1
    for motif_id, motif_seq in lib_seqs.items():
        d = levenshtein(var_region, motif_seq)
        if d < best_dist:
            best_dist = d
            best_id = motif_id
            if d == 0:
                break
    return best_id

def process_fastq(fastq_path, pool_name, lib_seqs):
    counts = defaultdict(Counter)
    n_total = n_amplicon = n_sample = n_motif = 0
    with gzip.open(fastq_path, "rt") as f:
        while True:
            header = f.readline()
            if not header: break
            seq = f.readline().strip().upper()
            f.readline(); f.readline()
            n_total += 1
            canonical, var = get_canonical_orientation(seq)
            if canonical is None: continue
            n_amplicon += 1
            sample = identify_sample(canonical, pool_name)
            if sample is None: continue
            n_sample += 1
            motif_id = assign_to_motif(var, lib_seqs)
            if motif_id is None: continue
            n_motif += 1
            counts[sample][motif_id] += 1
    print(f"\n{pool_name}:")
    print(f"  Total reads:               {n_total:,}")
    print(f"  Valid amplicon structure:  {n_amplicon:,}  ({n_amplicon/n_total*100:.1f}%)")
    print(f"  Assigned to a sample:      {n_sample:,}  ({n_sample/n_total*100:.1f}%)")
    print(f"  Assigned to a motif:       {n_motif:,}  ({n_motif/n_total*100:.1f}%)")
    for sample in counts:
        print(f"    {sample:12s} {sum(counts[sample].values()):,}")
    return counts


# ---------------------------------------------------------------------------
# STAGE 2: log2 fold-change enrichment
# ---------------------------------------------------------------------------

def compute_log2fc(all_counts, n_motifs=64):
    totals = {s: sum(c.values()) for s, c in all_counts.items()}
    def freq(s):
        return [(all_counts[s].get(i, 0) + PSEUDO) / (totals[s] + n_motifs * PSEUDO)
                for i in range(n_motifs)]
    lib_if = freq("Lib_input")
    dbd_if = freq("DBD_input")
    result = {}
    for i in range(n_motifs):
        lib_reps = [math.log2(freq(r)[i] / lib_if[i]) for r in ["Lib1","Lib2","Lib3"]]
        dbd_reps = [math.log2(freq(r)[i] / dbd_if[i]) for r in ["TF1","TF2","TF3"]]
        result[i] = {
            "lib_reps": lib_reps, "dbd_reps": dbd_reps,
            "lib_mean": float(np.mean(lib_reps)),
            "dbd_mean": float(np.mean(dbd_reps)),
            "delta": float(np.mean(dbd_reps) - np.mean(lib_reps)),
        }
    return result

def classify_tier(entry, motif_num):
    lib, dbd = entry["lib_reps"], entry["dbd_reps"]
    delta = entry["delta"]
    if min(lib) > 0.5 and min(dbd) > 0.5 and motif_num != 0:
        return "self_activator"
    if min(dbd) <= max(lib):
        return "none"
    if delta > T1_DELTA: return "T1"
    if delta > T2_DELTA: return "T2"
    return "none"


# ---------------------------------------------------------------------------
# STAGE 3: Gene annotation from BLAT + Table Browser
# ---------------------------------------------------------------------------

def load_blat_hits(path):
    with open(path) as f: return json.load(f)

def load_table_browser_tsv(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"): continue
            parts = line.split("\t")
            if len(parts) < 6: continue
            try:
                tx_start = int(float(parts[3]))
                tx_end   = int(float(parts[4]))
            except (ValueError, IndexError):
                continue
            rows.append({
                "refseq": parts[0].split(".")[0],
                "chrom": parts[1], "strand": parts[2],
                "tx_start": tx_start, "tx_end": tx_end,
                "gene": parts[5] if parts[5] else "UNKNOWN",
                "tss": tx_start if parts[2] == "+" else tx_end,
            })
    return rows


# ---------------------------------------------------------------------------
# STAGE 4: Odom 2004 tissue classification
# ---------------------------------------------------------------------------

def _norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s or "").upper())

def build_odom_lookup():
    hep_rs = {rs.split(".")[0] for _,rs in ODOM_S2_HEPATOCYTE}
    isl_rs = {rs.split(".")[0] for _,rs in ODOM_S3_ISLET}
    hep_nm = {_norm(g) for g,_ in ODOM_S2_HEPATOCYTE}
    isl_nm = {_norm(g) for g,_ in ODOM_S3_ISLET}
    return hep_rs, isl_rs, hep_nm, isl_nm

def classify_tissue(gene, refseq, odom):
    hep_rs, isl_rs, hep_nm, isl_nm = odom
    gn = _norm(gene)
    in_hep = (refseq in hep_rs) or (gn in hep_nm)
    in_isl = (refseq in isl_rs) or (gn in isl_nm)
    evidence = []
    if refseq in hep_rs: evidence.append("S2_refseq")
    elif gn in hep_nm: evidence.append("S2_name")
    if refseq in isl_rs: evidence.append("S3_refseq")
    elif gn in isl_nm: evidence.append("S3_name")
    if in_hep and in_isl: return "both", "|".join(evidence)
    if in_hep: return "liver", "|".join(evidence)
    if in_isl: return "islet", "|".join(evidence)
    return "neither", ""

def assign_gene(blat_hit, tb_rows, odom, window_bp=15000):
    if not blat_hit: return None
    chrom, pos = blat_hit["chrom"], blat_hit["start"]
    cands = []
    for r in tb_rows:
        if r["chrom"] != chrom: continue
        inside = r["tx_start"] <= pos <= r["tx_end"]
        tss_dist = abs(r["tss"] - pos)
        if inside or tss_dist <= window_bp:
            t, ev = classify_tissue(r["gene"], r["refseq"], odom)
            cands.append({"row": r, "inside": inside, "tss_dist": tss_dist,
                          "signed_dist": r["tss"] - pos, "tissue": t, "evidence": ev})
    if not cands:
        return {"gene": "intergenic", "tissue": "intergenic"}
    cands.sort(key=lambda c: c["tss_dist"])
    best = cands[0]
    return {
        "gene": best["row"]["gene"], "refseq": best["row"]["refseq"],
        "tss_dist": best["signed_dist"], "inside": best["inside"],
        "tissue": best["tissue"], "evidence": best["evidence"],
    }


# ---------------------------------------------------------------------------
# STAGE 5: Three-panel figure
# ---------------------------------------------------------------------------

TISSUE_COLOURS = {
    "liver": "#E85A4F", "islet": "#3B82BD", "both": "#8B4A9C",
    "neither": "#B8B8B8", "intergenic": "#E0E0E0", "unresolved": "#D0D0D0",
    "positive_control": "#2E8B57",
}
TIER_SHAPES = {"T1": ("*",320), "T2": ("D",130), "self_activator": ("X",150), "none": ("o",55)}

def make_3panel_plot(motifs, output_path):
    fig, axes = plt.subplots(3, 1, figsize=(19, 12), sharex=True)
    for ax, key, title, yl in [
        (axes[0], "lib_reps", "A. Lib-only (no HNF1A-DBD) - 3 replicates", "log$_2$FC vs input"),
        (axes[1], "dbd_reps", "B. DBD-Lib (HNF1A-DBD expressed) - 3 replicates", "log$_2$FC vs input"),
        (axes[2], None,       "C. $\\Delta$log$_2$FC (B$-$A) = HNF1A-DBD-specific enrichment", "$\\Delta$log$_2$FC"),
    ]:
        for i, m in enumerate(motifs):
            c = TISSUE_COLOURS.get(m["tissue"], "#888")
            sh, sz = TIER_SHAPES.get(m["tier"], ("o", 55))
            if key:
                ys = m[key]
                ax.scatter([i]*3, ys, c=c, marker=sh, s=sz*0.55, alpha=0.75,
                           edgecolors="k", linewidths=0.3)
                ax.hlines(np.mean(ys), i-0.32, i+0.32, colors=c, linewidth=2.2)
            else:
                ax.scatter([i], [m["delta"]], c=c, marker=sh, s=sz,
                           alpha=0.9, edgecolors="k", linewidths=0.5)
        ax.axhline(0, color="k", lw=0.5, ls="--", alpha=0.5)
        ax.set_ylabel(yl, fontsize=10)
        ax.set_title(title, fontsize=11, loc="left", pad=6)
        ax.grid(axis="y", alpha=0.2)
        ax.set_xlim(-1, len(motifs))

    axes[2].axhline(T1_DELTA, color="#E85A4F", lw=0.8, ls=":", alpha=0.7)
    axes[2].axhline(T2_DELTA, color="#D4A017", lw=0.8, ls=":", alpha=0.7)

    labels = []
    for m in motifs:
        g = m.get("gene") or ""
        if g in ("intergenic","ambiguous","HNF1_consensus_ctrl","?",""):
            labels.append("HNF1_ctrl" if m["motif_num"]==0 else "m"+str(m["motif_num"]))
        else:
            labels.append(g)
    axes[2].set_xticks(np.arange(len(motifs)))
    axes[2].set_xticklabels(labels, rotation=90, fontsize=6.5, ha="center")
    axes[2].set_xlabel("Nearest gene (ordered by DBD-Lib mean log$_2$FC, low to high)", fontsize=10)

    tc = Counter(m["tissue"] for m in motifs)
    leg1 = [
        mpatches.Patch(facecolor=TISSUE_COLOURS["liver"], ec="k",
                       label=f"Liver HNF1A target (Odom S2)  n={tc.get('liver',0)}"),
        mpatches.Patch(facecolor=TISSUE_COLOURS["islet"], ec="k",
                       label=f"Islet HNF1A target (Odom S3)  n={tc.get('islet',0)}"),
        mpatches.Patch(facecolor=TISSUE_COLOURS["both"], ec="k",
                       label=f"Both liver & islet  n={tc.get('both',0)}"),
        mpatches.Patch(facecolor=TISSUE_COLOURS["neither"], ec="k",
                       label=f"Not in Odom  n={tc.get('neither',0)}"),
        mpatches.Patch(facecolor=TISSUE_COLOURS["intergenic"], ec="k",
                       label=f"Intergenic/unresolved  n={tc.get('intergenic',0)+tc.get('unresolved',0)}"),
        mpatches.Patch(facecolor=TISSUE_COLOURS["positive_control"], ec="k",
                       label="Canonical HNF1 consensus (+ctrl)"),
    ]
    axes[0].legend(handles=leg1, loc="upper left", fontsize=8, framealpha=0.95, ncol=2)

    leg2 = [
        Line2D([0],[0], marker="*", color="w", mfc="#555", ms=14, mec="k",
               label=f"Tier 1 ($\\Delta$>{T1_DELTA}) n={sum(1 for m in motifs if m['tier']=='T1')}"),
        Line2D([0],[0], marker="D", color="w", mfc="#555", ms=9, mec="k",
               label=f"Tier 2 ({T2_DELTA}<$\\Delta$<{T1_DELTA}) n={sum(1 for m in motifs if m['tier']=='T2')}"),
        Line2D([0],[0], marker="X", color="w", mfc="#555", ms=10, mec="k",
               label=f"Self-activator n={sum(1 for m in motifs if m['tier']=='self_activator')}"),
        Line2D([0],[0], marker="o", color="w", mfc="#555", ms=6, mec="k",
               label=f"Non-hit n={sum(1 for m in motifs if m['tier']=='none')}"),
    ]
    axes[2].legend(handles=leg2, loc="upper left", fontsize=8, framealpha=0.95, ncol=2)

    plt.suptitle("HNF1A B1H Phase 1 - 63 promoter variants with Odom 2004 tissue annotation",
                 fontsize=13, y=0.995)
    plt.tight_layout()
    plt.subplots_adjust(top=0.94, hspace=0.28)
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved figure: {output_path}")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("HNF1A B1H Phase 1 - full pipeline")
    print("=" * 70)

    lib = load_library(INPUT_DIR / "hnf1a_63_unique_motif.txt")
    print(f"Loaded {len(lib)} library variants (including motif_00 positive control)")

    print("\nStage 1: demultiplex and map reads...")
    lib_counts = process_fastq(INPUT_DIR / "Lib-only_pool_Tube1_raw_fastq.gz", "Tube1", lib)
    dbd_counts = process_fastq(INPUT_DIR / "DBD-Lib_pool_Tube2_raw_fastq.gz", "Tube2", lib)
    all_counts = {**lib_counts, **dbd_counts}

    print("\nStage 2: enrichment calculations...")
    enrich = compute_log2fc(all_counts, n_motifs=64)
    for i, e in enrich.items():
        e["tier"] = classify_tier(e, i)
    n_t1 = sum(1 for e in enrich.values() if e["tier"] == "T1")
    n_t2 = sum(1 for e in enrich.values() if e["tier"] == "T2")
    print(f"  Tier 1 hits: {n_t1}")
    print(f"  Tier 2 hits: {n_t2}")

    print("\nStage 3: gene annotation from pre-computed BLAT + Table Browser...")
    blat = load_blat_hits(INPUT_DIR / "blat_best_hits.json")
    tb_rows = load_table_browser_tsv(INPUT_DIR / "BLAT_gene_overlaps.tsv")
    print(f"  BLAT resolved: {len(blat['resolved'])}")
    print(f"  BLAT ambiguous: {len(blat.get('ambiguous', {}))}")
    print(f"  Table Browser rows: {len(tb_rows)}")

    print("\nStage 4: Odom 2004 cross-reference...")
    odom = build_odom_lookup()
    print(f"  Odom S2 hepatocyte: {len(ODOM_S2_HEPATOCYTE)} entries")
    print(f"  Odom S3 islet:      {len(ODOM_S3_ISLET)} entries")

    motifs = []
    for i in range(64):
        if i == 8: continue  # library dropout, no reads
        mid = f"motif_{i:02d}"
        e = enrich[i]
        if i == 0:
            ann = {"gene": "HNF1_consensus_ctrl", "tissue": "positive_control",
                   "refseq": None, "tss_dist": None, "evidence": "Frain 1989"}
        else:
            blat_hit = blat["resolved"].get(mid)
            if blat_hit is None:
                ann = {"gene": "ambiguous", "tissue": "unresolved",
                       "refseq": None, "tss_dist": None, "evidence": ""}
            else:
                a = assign_gene(blat_hit, tb_rows, odom)
                ann = a if a else {"gene": "intergenic", "tissue": "intergenic"}
        motifs.append({
            "motif_id": mid, "motif_num": i,
            "gene": ann.get("gene", "?"), "refseq": ann.get("refseq"),
            "tissue": ann.get("tissue", "neither"),
            "tss_dist": ann.get("tss_dist"), "evidence": ann.get("evidence", ""),
            "lib_reps": e["lib_reps"], "dbd_reps": e["dbd_reps"],
            "lib_mean": e["lib_mean"], "dbd_mean": e["dbd_mean"],
            "delta": e["delta"], "tier": e["tier"],
        })

    motifs.sort(key=lambda m: m["dbd_mean"])

    # --- Manual annotation corrections (June 2026) ---
    # Two motifs needed manual curation:
    #
    # motif_37 (FCAMR, NM_001170631): The automated Odom lookup missed this
    # target because Odom 2004 listed it as FKSG87 (NM_032029), a legacy
    # placeholder for the same gene (NCBI Gene ID 83953). Both NM numbers
    # map to the same locus; manual verification confirmed S2 (hepatocyte)
    # classification. Reference: Odom et al., Science 2004, Table S2.
    #
    # motif_47 (ZNF44, NM_016264): BLAT returned two equal-score alignments
    # (chr19:12295831 and chr1:228850598), preventing automated gene
    # assignment. Manual inspection of the chr19 hit confirmed proximity to
    # ZNF44 (GIOT-2 in Odom 2004), present in both S2 and S3 tables.
    # Reference: Odom et al., Science 2004, Tables S2 and S3.
    MANUAL_CORRECTIONS = {
        37: {"gene": "FCAMR",  "refseq": "NM_001170631", "tissue": "liver",
             "odom_evidence": "S2_manual_FKSG87_alias"},
        47: {"gene": "ZNF44",  "refseq": "NM_016264",    "tissue": "both",
             "odom_evidence": "S2_S3_manual_GIOT2_alias"},
    }
    print("\nApplying manual annotation corrections:")
    for m in motifs:
        if m["motif_num"] in MANUAL_CORRECTIONS:
            corr = MANUAL_CORRECTIONS[m["motif_num"]]
            m.update(corr)
            print(f"  motif_{m['motif_num']:02d} -> gene={corr['gene']}, tissue={corr['tissue']}")
    # ---

    print("\nStage 5: figure + output tables...")
    make_3panel_plot(motifs, OUTPUT_DIR / "HNF1A_B1H_3panel_VERIFIED.png")

    with open(OUTPUT_DIR / "HNF1A_B1H_VERIFIED_ordered_with_hits.tsv","w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["rank","motif","gene","refseq","tissue","odom_evidence","tier",
                    "tss_dist_bp","lib_mean","dbd_mean","delta_log2fc"])
        for rank, m in enumerate(motifs, 1):
            w.writerow([rank, m["motif_id"], m["gene"] or "", m["refseq"] or "",
                        m["tissue"], m["evidence"], m["tier"],
                        m["tss_dist"] if m["tss_dist"] is not None else "",
                        f"{m['lib_mean']:.3f}", f"{m['dbd_mean']:.3f}",
                        f"{m['delta']:.3f}"])
    print(f"  Saved: {OUTPUT_DIR / 'HNF1A_B1H_VERIFIED_ordered_with_hits.tsv'}")

    tc = Counter(m["tissue"] for m in motifs)

    n_odom = sum(1 for m in motifs if m["tissue"] in ("liver","islet","both"))
    print(f"\n{n_odom}/{len(motifs)-1} library variants match Odom 2004 ({100*n_odom/(len(motifs)-1):.0f}%)")

    # Print Tier 1 and 2 hits
    print("\nTier 1 hits:")
    for m in sorted([m for m in motifs if m["tier"]=="T1"], key=lambda x: -x["delta"]):
        print(f"  {m['motif_id']}  {m['gene']:<14s} Δ={m['delta']:+.2f}  {m['tissue']}")
    print("\nTier 2 hits:")
    for m in sorted([m for m in motifs if m["tier"]=="T2"], key=lambda x: -x["delta"]):
        print(f"  {m['motif_id']}  {m['gene']:<14s} Δ={m['delta']:+.2f}  {m['tissue']}")

if __name__ == "__main__":
    main()
