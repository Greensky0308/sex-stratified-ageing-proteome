#!/usr/bin/env python3
"""Dump per-SNP z-scores of SomaScan CST3 and sarcopenia F/M across the cystatin locus."""
import os, sys
import pysam

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cst3_locus_inspect import read_somascan, read_sarc, SARC   # noqa: E402

LO, HI = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (23_540_000, 23_680_000)
ZMIN = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0


def main():
    ssc = read_somascan()
    sarc = {s: read_sarc(acc, ssc) for s, acc in SARC.items()}
    rows = []
    for rs, r in ssc.items():
        if not (LO <= r["pos"] <= HI):
            continue
        zc = r["beta"] / r["se"]
        out = [rs, r["pos"], zc]
        for s in ("F", "M"):
            out.append(sarc[s][rs]["beta"] / sarc[s][rs]["se"] if rs in sarc[s] else None)
        rows.append(out)
    rows.sort(key=lambda x: x[1])
    print(f"{'rsid':16s}{'pos':>10s}{'z_CST3':>10s}{'z_sarcF':>10s}{'z_sarcM':>10s}")
    for rs, pos, zc, zf, zm in rows:
        f = f"{zf:>10.1f}" if zf is not None else "         ."
        m = f"{zm:>10.1f}" if zm is not None else "         ."
        if abs(zc) > ZMIN or (zf is not None and abs(zf) > ZMIN):
            print(f"{rs:16s}{pos:>10d}{zc:>10.1f}{f}{m}")


if __name__ == "__main__":
    main()
