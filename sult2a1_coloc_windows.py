#!/usr/bin/env python3
"""SULT2A1 coloc across window sizes: SomaScan pQTL x sex-stratified UKB ALM.

Tests whether the default (+/-1 Mb) window inflates H3 by pulling in a distinct
ALM signal ~650 kb away; narrow windows centred on the gene should be used.
"""
import os, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coloc_abf import coloc                              # noqa: E402
from sult2a1_crossplatform import read_somascan, SOMA, CHR, LO, HI, COMP  # noqa: E402
from alm_locus_dump import dump                          # noqa: E402
from alm_sexdiff import ACC                              # noqa: E402

GENE = (47870327, 47886479)
MID = sum(GENE) // 2


def harm(b, ea_o, oa_o, ea_s, oa_s):
    if ea_s == ea_o and oa_s == oa_o:
        return b
    if ea_s == oa_o and oa_s == ea_o:
        return -b
    if ea_s == COMP.get(ea_o) and oa_s == COMP.get(oa_o):
        return b
    if ea_s == COMP.get(oa_o) and oa_s == COMP.get(ea_o):
        return -b
    return None


def main():
    ss = read_somascan(SOMA)
    alm = {s: dump(acc, CHR, LO, HI) for s, acc in ACC.items()}
    # harmonise ALM to SomaScan effect allele; keep only shared SNPs
    harm_alm = {}
    for s in ("F", "M"):
        h = {}
        for rs, (pos, ea_o, oa_o, b, se, z) in alm[s].items():
            if rs not in ss:
                continue
            bh = harm(b, ea_o, oa_o, ss[rs]["ea"], ss[rs]["oa"])
            if bh is None:
                continue
            h[rs] = dict(pos=pos, beta=bh, se=se)
        harm_alm[s] = h
        print(f"ALM {s}: harmonised SNPs in +/-1Mb = {len(h)}")

    print(f"\ngene={GENE}  MID={MID}")
    print(f"{'window':>10s}   {'ALM-F (n, H4)':>22s}   {'ALM-M (n, H4)':>22s}")
    for w in (25_000, 50_000, 100_000, 200_000, 500_000, 1_000_000):
        cells = []
        for s in ("F", "M"):
            b = {rs: (r["beta"], r["se"]) for rs, r in harm_alm[s].items()
                 if abs(r["pos"] - MID) <= w}
            a = {rs: (r["beta"], r["se"]) for rs, r in ss.items()
                 if abs(r["pos"] - MID) <= w and rs in b}
            n, pp = coloc(a, b)
            cells.append(f"n={n:5d} H3={pp[3]:.3f} H4={pp[4]:.3f}")
        print(f"{'+/-'+str(w//1000)+'kb':>10s}   {cells[0]:>22s}   {cells[1]:>22s}")


if __name__ == "__main__":
    main()
