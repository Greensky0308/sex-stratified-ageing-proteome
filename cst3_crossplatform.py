#!/usr/bin/env python3
"""CST3 cross-platform check (deCODE SomaScan pQTL -> sarcopenia F/M).

Steps
 1. Read deCODE SomaScan CST3 cis region (chr20 +/-1Mb of gene).
 2. Fetch the same region from the sex-stratified sarcopenia GWAS (bgzip+tabix)
    and harmonise betas onto the deCODE effect allele.
 3. Wald-ratio MR per sex for each cis instrument; sex-difference z.
 4. Write coloc.abf input tables (rsid pos beta se).

deCODE cols: Chrom Pos Name rsids effectAllele otherAllele Beta Pval minus_log10_pval SE N ImpMAF
sarcopenia cols: chromosome base_pair_location effect_allele other_allele beta standard_error
                 effect_allele_frequency p_value rsid ...
"""
import math, os, sys
import pysam

ROOT = os.path.dirname(os.path.abspath(__file__))
CIS_RAW = os.environ.get("CST3_CIS_RAW") or os.path.join(ROOT, "out", "cst3_cis_raw.tsv")
SARC = {"F": "GCST90832979", "M": "GCST90832980"}
REGION = ("20", 22626706, 24638556)
OUT = os.path.join(ROOT, "out")
INSTRUMENTS = ["rs911119", "rs6114209"]          # UKB-PPP lead + SomaScan lead
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def read_somascan(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 12:
                continue
            rs = c[3].split(",")[0].strip()
            if not rs.startswith("rs"):
                continue
            try:
                d[rs] = dict(pos=int(c[1]), ea=c[4], oa=c[5], beta=float(c[6]),
                             se=float(c[9]), p=float(c[7]))
            except ValueError:
                continue
    return d


def harmonise(beta, ea_o, oa_o, ea_s, oa_s):
    """Return beta_o expressed per SomaScan effect allele, or None if incompatible."""
    if ea_s == ea_o and oa_s == oa_o:
        return beta
    if ea_s == oa_o and oa_s == ea_o:
        return -beta
    # strand flip
    ea_c = COMP.get(ea_o); oa_c = COMP.get(oa_o)
    if ea_s == ea_c and oa_s == oa_c:
        return beta
    if ea_s == oa_c and oa_s == ea_c:
        return -beta
    return None


def fetch_outcome(acc, ssc):
    path = os.path.join(ROOT, "data", "sumstats", f"{acc}.h.tsv.gz")
    tb = pysam.TabixFile(path)
    rows = {}
    crm, lo, hi = REGION
    n_tot = n_ok = 0
    for row in tb.fetch(crm, lo, hi):
        c = row.split("\t")
        if len(c) < 9:
            continue
        rs = c[8]
        n_tot += 1
        if rs not in ssc:
            continue
        s = ssc[rs]
        try:
            b = float(c[4]); se = float(c[5])
        except ValueError:
            continue
        bh = harmonise(b, c[2], c[3], s["ea"], s["oa"])
        if bh is None:
            continue
        rows[rs] = dict(pos=int(c[1]), beta=bh, se=se, ea=s["ea"])
        n_ok += 1
    tb.close()
    print(f"  {acc}: region rows={n_tot}, harmonised to SomaScan={n_ok}")
    return rows


def wald(bx, sx, by, sy):
    if abs(bx) < 1e-9 or sy <= 0:
        return None
    return by / bx, sy / abs(bx)


def main():
    os.makedirs(OUT, exist_ok=True)
    ssc = read_somascan(CIS_RAW)
    print(f"SomaScan cis SNPs: {len(ssc)}")

    # MR per instrument
    print("\n=== Wald-ratio MR (SomaScan instrument -> sarcopenia) ===")
    res = {}
    for sex, acc in SARC.items():
        out = fetch_outcome(acc, ssc)
        res[sex] = out
        # write coloc input
        with open(os.path.join(OUT, f"cst3_coloc_sarc{sex}.tsv"), "w") as fh:
            for rs in sorted(out):
                r = out[rs]
                fh.write(f"{rs}\t{r['pos']}\t{r['beta']}\t{r['se']}\n")

    with open(os.path.join(OUT, "cst3_coloc_somascan.tsv"), "w") as fh:
        for rs in sorted(ssc):
            r = ssc[rs]
            fh.write(f"{rs}\t{r['pos']}\t{r['beta']}\t{r['se']}\n")

    hdr = "instrument\trs\tsomascan_beta\tbF\tseF\tbM\tseM\tdiff\tse_diff\tzdiff\tp"
    print("\n" + hdr)
    for rs in INSTRUMENTS:
        if rs not in ssc:
            print(f"{rs}: NOT in SomaScan cis"); continue
        bx = ssc[rs]["beta"]; sx = ssc[rs]["se"]
        vals = {}
        for sex in ("F", "M"):
            if rs not in res[sex]:
                vals[sex] = None; continue
            o = res[sex][rs]
            vals[sex] = wald(bx, sx, o["beta"], o["se"])
        if vals["F"] and vals["M"]:
            (bF, seF), (bM, seM) = vals["F"], vals["M"]
            diff = bF - bM
            sd = math.sqrt(seF**2 + seM**2)
            z = diff / sd
            p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
            print(f"CST3\t{rs}\t{bx:+.4f}\t{bF:+.3f}\t{seF:.3f}\t{bM:+.3f}\t{seM:.3f}\t"
                  f"{diff:+.3f}\t{sd:.3f}\t{z:+.3f}\t{p:.3g}")

    # cross-check against the UKB-PPP numbers
    print("\n[ref UKB-PPP rs911119] bF=-14.775 seF=0.154 bM=-15.821 seM=0.173 "
          "diff=+1.045 zdiff=+4.516 p=6.3e-6 q=0.011")
    print("(sex-diff z is invariant to exposure scaling when the SAME SNP is used)")


if __name__ == "__main__":
    main()
