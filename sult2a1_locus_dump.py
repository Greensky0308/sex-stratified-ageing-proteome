#!/usr/bin/env python3
"""Side-by-side z-scores of SomaScan SULT2A1 pQTL and ALM F/M across the locus."""
import os, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sult2a1_crossplatform import read_somascan, SOMA, CHR, LO, HI, COMP   # noqa: E402
from alm_locus_dump import dump                                           # noqa: E402
from alm_sexdiff import ACC                                               # noqa: E402

ZMIN = float(sys.argv[1]) if len(sys.argv) > 1 else 4.0


def main():
    ss = read_somascan(SOMA)
    F = dump(ACC["F"], CHR, LO, HI)
    M = dump(ACC["M"], CHR, LO, HI)

    def orient(raw, rs):
        pos, ea_o, oa_o, b, se, z = raw[rs]
        ea_s = ss[rs]["ea"]
        if ea_o == ea_s:
            return z
        if oa_o == ea_s:
            return -z
        return None

    rows = []
    for rs in ss:
        if rs not in F or rs not in M:
            continue
        zs = ss[rs]["beta"] / ss[rs]["se"]
        zf = orient(F, rs)
        zm = orient(M, rs)
        if zf is None or zm is None:
            continue
        rows.append((ss[rs]["pos"], rs, zs, zf, zm, ss[rs]["p"]))
    rows.sort()
    print(f"{'rsid':16s}{'pos':>10s}{'z_pQTL':>9s}{'z_ALMF':>9s}{'z_ALMM':>9s}")
    for pos, rs, zs, zf, zm, p in rows:
        if abs(zs) > ZMIN or abs(zf) > ZMIN or abs(zm) > ZMIN:
            print(f"{rs:16s}{pos:>10d}{zs:>9.1f}{zf:>9.1f}{zm:>9.1f}")

    import numpy as np
    zs = np.array([r[2] for r in rows]); zf = np.array([r[3] for r in rows]); zm = np.array([r[4] for r in rows])
    print(f"\nn={len(rows)}  corr(z_pQTL, z_ALMF)={np.corrcoef(zs,zf)[0,1]:.3f}  "
          f"corr(z_pQTL, z_ALMM)={np.corrcoef(zs,zm)[0,1]:.3f}")
    for lbl, z in (("F", zf), ("M", zm)):
        top = max(range(len(rows)), key=lambda i: abs(z[i]))
        print(f"ALM {lbl} peak: {rows[top][1]} pos={rows[top][0]} z={z[top]:+.1f}")
    top = max(range(len(rows)), key=lambda i: abs(zs[i]))
    print(f"pQTL peak:   {rows[top][1]} pos={rows[top][0]} z={zs[top]:+.1f}")


if __name__ == "__main__":
    main()
