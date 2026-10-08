#!/usr/bin/env python3
"""Sex-differential cis-pQTL MR on Neale-lab UKB whole-body fat-free mass (field 23104).

Neale round2 files are GRCh37, no rsid, columns:
  variant(chr:pos:ref:alt) minor_allele minor_AF low_confidence_variant n_complete_samples
  AC ytx beta se tstat pval
Instruments (data/instruments_ukbppp_grch37.tsv) carry GRCh37 coordinates too -> match on chr:pos.

Allele orientation is NOT resolved (no stable effect-allele column); the F and M files
share one convention, so the sex-difference test is orientation-invariant (a global
flip cancels). Absolute sign is reported relative to that convention.
"""
import argparse, math, os
import numpy as np
from scipy.stats import norm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_inst(path):
    rows = []
    with open(path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            try:
                rows.append((c[0], c[1], str(c[2]), int(c[3]), c[4], float(c[5])))
            except ValueError:
                continue
    return rows


def extract(inst, path):
    key = {(ch, pos): (prot, rs, ea, bx) for prot, rs, ch, pos, ea, bx in inst if abs(bx) > 1e-9}
    import gzip
    out = {}
    with gzip.open(path, "rt") as fh:
        cm = fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            v = c[0].split(":")
            if len(v) < 4:
                continue
            try:
                t = key.get((v[0], int(v[1])))
            except ValueError:
                continue
            if t is None:
                continue
            prot, rs, ea, bx = t
            try:
                b = float(c[7]); se = float(c[8])
            except ValueError:
                continue
            if se <= 0:
                continue
            out[prot] = (b / bx, se / abs(bx))
    return out


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[o[i]] * n / (i + 1)); q[o[i]] = prev
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--field", default="23101", help="Neale field id, e.g. 23101 whole-body FFM")
    ap.add_argument("--instruments", default=os.path.join(ROOT, "data", "instruments_ukbppp_grch37.tsv"))
    a = ap.parse_args()
    files = {"F": f"neale_{a.field}_female.tsv.bgz", "M": f"neale_{a.field}_male.tsv.bgz"}
    inst = load_inst(a.instruments)
    per = {}
    for sex, fn in files.items():
        path = os.path.join(ROOT, "data", "sumstats", fn)
        per[sex] = extract(inst, path)
        print(f"  {fn}: matched {len(per[sex])}")
    res = []
    for prot in per["F"]:
        if prot not in per["M"]:
            continue
        bF, sF = per["F"][prot]; bM, sM = per["M"][prot]
        if sF <= 0 or sM <= 0:
            continue
        res.append((prot, bF, sF, bM, sM, (bF - bM) / math.sqrt(sF**2 + sM**2)))
    p = [2 * (1 - norm.cdf(abs(r[5]))) for r in res]
    q = bh(p)
    order = sorted(range(len(res)), key=lambda i: abs(res[i][5]), reverse=True)
    outdir = os.path.join(ROOT, "out", f"neale_{a.field}")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "ffm_sexdiff.tsv"), "w") as fh:
        fh.write("protein\tbF\tseF\tbM\tseM\tzdiff\tp\tq\n")
        for i in order:
            r = res[i]
            fh.write(f"{r[0]}\t{r[1]:.4g}\t{r[2]:.4g}\t{r[3]:.4g}\t{r[4]:.4g}\t{r[5]:.3f}\t{p[i]:.3g}\t{q[i]:.3g}\n")
    print(f"paired={len(res)}  q<0.05={int((np.array(q)<0.05).sum())}")
    for i in order[:10]:
        r = res[i]
        print(f"  {r[0]:12s} bF={r[1]:+.4f} bM={r[3]:+.4f} z={r[5]:+.2f} p={p[i]:.2e} q={q[i]:.3g}")
    for i, r in enumerate(res):
        if r[0] == "SULT2A1":
            print(f"[SULT2A1] bF={r[1]:+.4f} bM={r[3]:+.4f} z={r[5]:+.2f} p={p[i]:.3g} q={q[i]:.3g}")


if __name__ == "__main__":
    main()
