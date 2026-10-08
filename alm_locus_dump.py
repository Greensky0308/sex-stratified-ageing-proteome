#!/usr/bin/env python3
"""Dump per-SNP z-scores of ALM female/male in a genomic window (to inspect a locus).

Usage: alm_locus_dump.py <chrom> <start> <end> [zmin]
"""
import gzip, os, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from alm_sexdiff import _colmap, ACC   # noqa: E402


def dump(acc, chrom, lo, hi):
    path = os.path.join(ROOT, "data", "sumstats", f"{acc}.h.tsv.gz")
    rows = {}
    with gzip.open(path, "rt") as fh:
        cm = _colmap(fh.readline())
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if c[cm["chrom"]] != chrom:
                continue
            try:
                pos = int(c[cm["pos"]])
            except ValueError:
                continue
            if pos < lo or pos > hi:
                continue
            try:
                b = float(c[cm["beta"]]); se = float(c[cm["se"]])
            except ValueError:
                continue
            if se <= 0:
                continue
            rows[c[cm["rsid"]]] = (pos, c[cm["ea"]], c[cm["oa"]], b, se, b / se)
    return rows


def main():
    chrom = sys.argv[1]; lo = int(sys.argv[2]); hi = int(sys.argv[3])
    zmin = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0
    F = dump(ACC["F"], chrom, lo, hi)
    M = dump(ACC["M"], chrom, lo, hi)
    print(f"window {chrom}:{lo}-{hi}  nF={len(F)} nM={len(M)}")
    print("zM is oriented onto the female effect allele; sign mismatch = real direction difference")
    print(f"{'rsid':16s}{'pos':>10s}{'ea':>4s}{'zF':>8s}{'zM':>8s}")
    rows = []
    for rs in set(F) | set(M):
        pos = F[rs][0] if rs in F else M[rs][0]
        zf = F[rs][5] if rs in F else None
        zm = None
        if rs in M and rs in F:
            mea, moa, mb, mse = M[rs][1], M[rs][2], M[rs][3], M[rs][4]
            fea = F[rs][1]
            if mea == fea:
                ori = mb
            elif moa == fea:
                ori = -mb
            else:
                ori = None
            if ori is not None:
                zm = ori / mse
        elif rs in M:
            zm = M[rs][5]
        rows.append((pos, rs, F[rs][1] if rs in F else "?", zf, zm))
    rows.sort()
    for pos, rs, ea, zf, zm in rows:
        if (zf is not None and abs(zf) > zmin) or (zm is not None and abs(zm) > zmin) \
           or rs == "rs62129966":
            a = f"{zf:>8.1f}" if zf is not None else "       ."
            b = f"{zm:>8.1f}" if zm is not None else "       ."
            print(f"{rs:16s}{pos:>10d}{ea:>4s}{a}{b}")


if __name__ == "__main__":
    main()
