#!/usr/bin/env python3
"""Sex-differential cis-pQTL MR on real appendicular lean mass (UKB, Pei 2020).

Replaces the circular "sarcopenia index" outcome (SI = creatinine/cystatin C).
Instruments: UKB-PPP top cis pQTL per protein (all cis by construction).
Outcomes: GCST90000027 (female ALM), GCST90000026 (male ALM).

Wald ratio per sex -> sex-difference z -> BH FDR.
"""
import argparse, math, os
import numpy as np
from scipy.stats import norm

ACC = {"F": "GCST90000027", "M": "GCST90000026"}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def load_instruments(path):
    rows = []
    with open(path) as fh:
        hdr = fh.readline().split("\t")
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            try:
                rows.append((c[0], c[1], str(c[2]), int(c[3]), c[4], float(c[5])))
            except ValueError:
                continue
    return rows


def _wald(prot, rsid, ea, bx, c):
    """c = outcome row [chrom,pos,ea,oa,beta,se,eaf,p,rsid,...]"""
    o_ea, o_oa = c[2], c[3]
    if ea == o_ea:
        sign = 1
    elif ea == o_oa:
        sign = -1
    elif ea == COMP.get(o_ea):
        sign = 1
    elif ea == COMP.get(o_oa):
        sign = -1
    else:
        return None
    return (prot, rsid, sign * float(c[4]) / bx, float(c[5]) / abs(bx))


def _colmap(header):
    h = [x.strip() for x in header.split("\t")]

    def idx(*names):
        for n in names:
            if n in h:
                return h.index(n)
        return None
    return {"rsid": idx("rsid", "hm_rsid", "variant_id", "hm_variant_id"),
            "ea": idx("effect_allele", "hm_effect_allele"),
            "oa": idx("other_allele", "hm_other_allele"),
            "beta": idx("beta", "hm_beta"),
            "se": idx("standard_error"),
            "chrom": idx("chromosome", "hm_chrom"),
            "pos": idx("base_pair_location", "hm_pos")}


def extract(inst, acc, root):
    """Stream-scan by rsid (instrument coordinates are GRCh37, outcomes GRCh38)."""
    import gzip
    path = os.path.join(root, "data", "sumstats", f"{acc}.h.tsv.gz")
    by_rs = {rsid: (prot, rsid, ea, bx)
             for (prot, rsid, ch, pos, ea, bx) in inst if abs(bx) > 1e-9}
    out, nrow, nmatched = [], 0, 0
    with gzip.open(path, "rt") as fh:
        cm = _colmap(fh.readline())
        for line in fh:
            nrow += 1
            c = line.rstrip("\n").split("\t")
            if len(c) <= max(v for v in cm.values() if v is not None):
                continue
            t = by_rs.get(c[cm["rsid"]])
            if t is None:
                continue
            nmatched += 1
            prot, rsid, ea, bx = t
            w = _wald(prot, rsid, ea, bx, [None, None, c[cm["ea"]], c[cm["oa"]],
                                           c[cm["beta"]], c[cm["se"]]])
            if w:
                out.append(w)
    print(f"  [{acc}] rows={nrow} matched={nmatched} wald={len(out)} cols={cm}")
    return out


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[o[i]] * n / (i + 1)); q[o[i]] = prev
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--instruments",
                    default=os.path.join(os.path.dirname(os.path.dirname(
                        os.path.abspath(__file__))), "data", "instruments_ukbppp_grch37.tsv"))
    ap.add_argument("--outdir", default="out/alm_ukbppp")
    ap.add_argument("--f-acc", default=ACC["F"])
    ap.add_argument("--m-acc", default=ACC["M"])
    ap.add_argument("--max-z", type=float, default=None)
    a = ap.parse_args()
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.makedirs(os.path.join(root, a.outdir), exist_ok=True)
    inst = load_instruments(a.instruments)
    print("instruments:", len(inst))

    per = {"F": {}, "M": {}}
    for sex, acc in {"F": a.f_acc, "M": a.m_acc}.items():
        rows = extract(inst, acc, root)
        print(f"  {acc} ({sex}): {len(rows)} MR rows")
        for prot, rsid, b, se in rows:
            per[sex][prot] = (b, se)

    res = []
    for prot in per["F"]:
        if prot not in per["M"]:
            continue
        bF, sF = per["F"][prot]; bM, sM = per["M"][prot]
        if sF <= 0 or sM <= 0:
            continue
        z = (bF - bM) / math.sqrt(sF**2 + sM**2)
        res.append((prot, bF, sF, bM, sM, z))
    if not res:
        print("no paired proteins"); return
    p = [2 * (1 - norm.cdf(abs(r[5]))) for r in res]
    q = bh(p)
    order = sorted(range(len(res)), key=lambda i: abs(res[i][5]), reverse=True)
    out = os.path.join(root, a.outdir, "alm_sexdiff.tsv")
    with open(out, "w") as fh:
        fh.write("protein\tbF\tseF\tbM\tseM\tzdiff\tp\tq\n")
        for i in order:
            r = res[i]
            fh.write(f"{r[0]}\t{r[1]:.4g}\t{r[2]:.4g}\t{r[3]:.4g}\t{r[4]:.4g}\t{r[5]:.3f}\t{p[i]:.3g}\t{q[i]:.3g}\n")
    n_sig = int((np.array(q) < 0.05).sum())
    print(f"paired proteins: {len(res)} | q<0.05: {n_sig}")
    print("top 12 by |z|:")
    for i in order[:12]:
        r = res[i]
        print(f"  {r[0]:12s} bF={r[1]:+.4f} bM={r[3]:+.4f} z={r[5]:+.2f} p={p[i]:.2e} q={q[i]:.3g}")
    print("wrote", out)
    for prot in ("CST3",):
        for i, r in enumerate(res):
            if r[0] == prot:
                print(f"[{prot}] bF={r[1]:+.4f} bM={r[3]:+.4f} z={r[5]:+.2f} p={p[i]:.3g} q={q[i]:.3g}")


if __name__ == "__main__":
    main()
