#!/usr/bin/env python3
"""Inspect the CST3 cis locus: lead SNPs on each side + coloc at multiple window sizes.

Reads the SomaScan CST3 cis extraction and the sarcopenia F/M region, then asks
whether H3=1 is driven by an over-wide window.
"""
import os, sys
import pysam

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coloc_abf import coloc, load                       # noqa: E402

CIS_RAW = os.environ.get("CST3_CIS_RAW") or os.path.join(ROOT, "out", "cst3_cis_raw.tsv")
SARC = {"F": "GCST90832979", "M": "GCST90832980"}
GENE = (23626706, 23638556)
MID = sum(GENE) // 2
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def read_somascan():
    d = {}
    for line in open(CIS_RAW):
        c = line.rstrip("\n").split("\t")
        if len(c) < 12:
            continue
        rs = c[3].split(",")[0].strip()
        if not rs.startswith("rs"):
            continue
        try:
            d[rs] = dict(pos=int(c[1]), beta=float(c[6]), se=float(c[9]), p=float(c[7]),
                         ea=c[4], oa=c[5])
        except ValueError:
            continue
    return d


def harm(b, ea_o, oa_o, ea_s, oa_s):
    if ea_s == ea_o and oa_s == oa_o: return b
    if ea_s == oa_o and oa_s == ea_o: return -b
    ec, oc = COMP.get(ea_o), COMP.get(oa_o)
    if ea_s == ec and oa_s == oc: return b
    if ea_s == oc and oa_s == ec: return -b
    return None


def read_sarc(acc, ssc):
    tb = pysam.TabixFile(os.path.join(ROOT, "data", "sumstats", f"{acc}.h.tsv.gz"))
    d = {}
    for row in tb.fetch("20", MID - 1_500_000, MID + 1_500_000):
        c = row.split("\t")
        if len(c) < 9 or c[8] not in ssc:
            continue
        s = ssc[c[8]]
        try:
            b = float(c[4]); se = float(c[5]); p = float(c[7])
        except ValueError:
            continue
        bh = harm(b, c[2], c[3], s["ea"], s["oa"])
        if bh is None:
            continue
        d[c[8]] = dict(pos=int(c[1]), beta=bh, se=se, p=p)
    tb.close()
    return d


def top(d, k=5):
    return sorted(d.items(), key=lambda kv: kv[1]["p"])[:k]


def main():
    ssc = read_somascan()
    print("=== CST3 (SomaScan) top 5 ===")
    for rs, r in top(ssc):
        print(f"  {rs:16s} pos={r['pos']} beta={r['beta']:+.4f} p={r['p']:.2e}")
    sarc = {s: read_sarc(acc, ssc) for s, acc in SARC.items()}
    for s in ("F", "M"):
        print(f"=== sarcopenia {s} top 5 in locus ===")
        for rs, r in top(sarc[s]):
            print(f"  {rs:16s} pos={r['pos']} beta={r['beta']:+.4f} p={r['p']:.2e}")

    a = {rs: (r["beta"], r["se"]) for rs, r in ssc.items()}
    print("\n=== coloc at varying window (SomaScan CST3 x sarcopenia) ===")
    for w in (50_000, 100_000, 250_000, 500_000, 1_000_000):
        line = f"  +/-{w//1000:>4d}kb :"
        for s in ("F", "M"):
            b = {rs: (r["beta"], r["se"]) for rs, r in sarc[s].items()
                 if abs(r["pos"] - MID) <= w}
            aa = {rs: v for rs, v in a.items() if abs(ssc[rs]["pos"] - MID) <= w}
            n, pp = coloc(aa, b)
            line += f"  [sarc{s} n={n} H3={pp[3]:.3f} H4={pp[4]:.3f}]"
        print(line)


if __name__ == "__main__":
    main()
