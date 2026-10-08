#!/usr/bin/env python3
"""SULT2A1 cross-platform check: deCODE SomaScan pQTL vs sex-stratified UKB ALM.

 - cis region of SULT2A1 (chr19 +/-1Mb) from the deCODE SomaScan file
 - ALM region from GCST90000027 (F) and GCST90000026 (M), harmonised to the
   SomaScan effect allele
 - Wald-ratio MR per sex (SomaScan instrument) + sex-difference z
 - coloc.abf SomaScan x ALM (F, M)
"""
import gzip, math, os, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coloc_abf import coloc                      # noqa: E402
from alm_locus_dump import dump                  # noqa: E402
from alm_sexdiff import ACC                      # noqa: E402

SOMA = os.path.join(ROOT, "data", "decode",
                    "Proteomics_SMP_PC0_9829_91_SULT2A1_SULT_2A1_10032022.txt.gz")
GENE = (47870327, 47886479)
CHR, LO, HI = "19", 46870327, 48886479
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def read_somascan(path):
    d = {}
    with gzip.open(path, "rt") as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 12 or c[0] != f"chr{CHR}":
                continue
            try:
                pos = int(c[1])
            except ValueError:
                continue
            if pos < LO or pos > HI:
                continue
            rs = c[3].split(",")[0].strip()
            if not rs.startswith("rs"):
                continue
            try:
                d[rs] = dict(pos=pos, ea=c[4], oa=c[5], beta=float(c[6]),
                             se=float(c[9]), p=float(c[7]))
            except ValueError:
                continue
    return d


def harm(beta, ea_o, oa_o, ea_s, oa_s):
    if ea_s == ea_o and oa_s == oa_o: return beta
    if ea_s == oa_o and oa_s == ea_o: return -beta
    if ea_s == COMP.get(ea_o) and oa_s == COMP.get(oa_o): return beta
    if ea_s == COMP.get(oa_o) and oa_s == COMP.get(ea_o): return -beta
    return None


def main():
    path = SOMA
    ss = read_somascan(path)
    print(f"SomaScan cis SNPs: {len(ss)}  (from {os.path.basename(path)})")
    lead_ss = min(ss, key=lambda r: ss[r]["p"])
    print(f"SomaScan lead: {lead_ss} pos={ss[lead_ss]['pos']} beta={ss[lead_ss]['beta']:+.4f} p={ss[lead_ss]['p']:.2e}")

    alm = {}
    for sex, acc in ACC.items():
        raw = dump(acc, CHR, LO, HI)
        h = {}
        for rs, (pos, ea_o, oa_o, b, se, z) in raw.items():
            if rs not in ss:
                continue
            bh = harm(b, ea_o, oa_o, ss[rs]["ea"], ss[rs]["oa"])
            if bh is None:
                continue
            h[rs] = dict(pos=pos, beta=bh, se=se)
        alm[sex] = h
        print(f"  ALM {sex}: region SNPs harmonised to SomaScan = {len(h)}")

    # MR per instrument
    print("\n=== Wald-ratio MR (SomaScan SULT2A1 -> ALM) ===")
    print(f"{'instrument':16s}{'bx':>9s}{'bF':>9s}{'bM':>9s}{'zdiff':>8s}{'p':>10s}")
    instruments = [lead_ss, "rs62129966"]
    for rs in instruments:
        if rs not in ss:
            print(f"{rs}: not in SomaScan"); continue
        bx = ss[rs]["beta"]
        v = {}
        for sex in ("F", "M"):
            v[sex] = alm[sex][rs]["beta"] / bx if rs in alm[sex] else None
        if v["F"] and v["M"]:
            z = (v["F"] - v["M"]) / math.sqrt((alm["F"][rs]["se"] / abs(bx))**2 +
                                              (alm["M"][rs]["se"] / abs(bx))**2)
            p = math.erfc(abs(z) / math.sqrt(2))
            print(f"{rs:16s}{bx:>9.4f}{v['F']:>9.4f}{v['M']:>9.4f}{z:>8.2f}{p:>10.2e}")

    # coloc
    ss_tab = {rs: (r["beta"], r["se"]) for rs, r in ss.items()}
    print("\n=== coloc.abf (SomaScan SULT2A1 x ALM) ===")
    for sex in ("F", "M"):
        t = {rs: (r["beta"], r["se"]) for rs, r in alm[sex].items()}
        n, pp = coloc(ss_tab, t)
        print(f"  ALM {sex}: nSNP={n} H3={pp[3]:.3f} H4(shared)={pp[4]:.3f}")


if __name__ == "__main__":
    main()
