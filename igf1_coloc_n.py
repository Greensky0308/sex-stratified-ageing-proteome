#!/usr/bin/env python3
"""Report the number of SNPs entering each IGF-1 colocalization window.

Reproduces the SNP count used by igf1_coloc.py, reading the local SomaScan
SULT2A1 file and the local Neale IGF-1 (field 30770) files instead of the
remote copies.
"""
import gzip, os, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
from sult2a1_crossplatform import read_somascan, SOMA, COMP   # noqa: E402
from pyliftover import LiftOver                                # noqa: E402

LO37, HI37 = 47_600_000, 48_900_000
MID38 = 47_871_693
WINDOWS = (1_300_000, 500_000, 250_000, 100_000, 50_000)


def read_igf1_region(path):
    """chr19 region from a Neale round-2 file: variant = chr:pos:ref:alt."""
    out = {}
    with gzip.open(path, "rt") as fh:
        fh.readline()
        for line in fh:
            if not line.startswith("19:"):
                continue
            c = line.split("\t")
            mk = c[0].split(":")
            p = int(mk[1])
            if p > HI37:
                break
            if p < LO37:
                continue
            try:
                out[p] = (mk[2], mk[3], float(c[7]), float(c[8]))
            except (ValueError, IndexError):
                continue
    return out


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
    bypos = {}
    for rs, r in ss.items():
        bypos.setdefault(r["pos"], []).append(r)
    lo = LiftOver("hg19", "hg38")
    for sex in ("female", "male"):
        reg = read_igf1_region(os.path.join(ROOT, "data", "igf1", f"igf1_{sex}.tsv.bgz"))
        keep = {}
        for p37, (ref, alt, b, se) in reg.items():
            res = lo.convert_coordinate("chr19", p37)
            if not res:
                continue
            p38 = res[0][1]
            if p38 not in bypos:
                continue
            for r in bypos[p38]:
                if harm(b, ref, alt, r["ea"], r["oa"]) is not None:
                    keep[p38] = True
                    break
        print(f"IGF-1 {sex}: region SNPs(GRCh37)={len(reg)}  harmonised={len(keep)}")
        for w in WINDOWS:
            n = sum(1 for p in keep if abs(p - MID38) <= w)
            print(f"   +/-{w//1000:>4d}kb  n={n}")


if __name__ == "__main__":
    main()
