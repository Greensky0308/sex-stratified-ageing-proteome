#!/usr/bin/env python3
"""Dump muscle-weakness (grip) F/M z-scores at the SULT2A1 locus, allele-oriented."""
import os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from alm_locus_dump import dump      # noqa: E402

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
F_ACC, M_ACC = "GCST90007527", "GCST90007528"
CHR = "19"


def main():
    lo, hi = int(sys.argv[1]), int(sys.argv[2])
    zmin = float(sys.argv[3]) if len(sys.argv) > 3 else 2.5
    F = dump(F_ACC, CHR, lo, hi)
    M = dump(M_ACC, CHR, lo, hi)
    print(f"chr{CHR}:{lo}-{hi}  nF={len(F)} nM={len(M)}  (zM oriented to female effect allele)")
    print(f"{'rsid':16s}{'pos':>10s}{'ea':>4s}{'zF':>8s}{'zM':>8s}")
    rows = []
    for rs in set(F) | set(M):
        pos = F[rs][0] if rs in F else M[rs][0]
        zf = F[rs][5] if rs in F else None
        zm = None
        if rs in F and rs in M:
            mea, moa, mb, mse = M[rs][1], M[rs][2], M[rs][3], M[rs][4]
            fea = F[rs][1]
            ori = mb if mea == fea else (-mb if moa == fea else None)
            if ori is not None:
                zm = ori / mse
        elif rs in M:
            zm = M[rs][5]
        rows.append((pos, rs, F[rs][1] if rs in F else "?", zf, zm))
    rows.sort()
    for pos, rs, ea, zf, zm in rows:
        if rs == "rs62129966" or (zf is not None and abs(zf) > zmin) or (zm is not None and abs(zm) > zmin):
            a = f"{zf:>8.2f}" if zf is not None else "       ."
            b = f"{zm:>8.2f}" if zm is not None else "       ."
            print(f"{rs:16s}{pos:>10d}{ea:>4s}{a}{b}")


if __name__ == "__main__":
    main()
