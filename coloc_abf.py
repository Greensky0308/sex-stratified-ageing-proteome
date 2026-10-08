#!/usr/bin/env python3
"""coloc.abf (Giambartolomei 2014) — pure-Python implementation, log-space.

Usage: coloc_abf.py <trait1.tsv> <trait2.tsv>   (cols: rsid pos beta se  [maf])
Priors p1=p2=1e-4, p12=1e-5, sd.prior(W)=0.15 (quantitative trait).

Works in log space so very strong signals (z ~ 1e2) do not overflow exp().
"""
import sys, math

NEG = float("-inf")


def load(p):
    d = {}
    for i, l in enumerate(open(p)):
        c = l.rstrip("\n").split("\t")
        if i == 0 and not c[1].replace('.', '').replace('-', '').isdigit():
            continue
        try:
            d[c[0]] = (float(c[2]), float(c[3]))
        except Exception:
            pass
    return d


W = 0.15 ** 2


def labf(b, se):
    V = se ** 2
    if V <= 0:
        return NEG
    r = W / (W + V)
    z = b / se
    return 0.5 * (math.log(1 - r) + r * z * z)


def logsumexp(xs):
    m = max(xs)
    if m == NEG:
        return NEG
    return m + math.log(sum(math.exp(x - m) for x in xs))


def log_diff(a, b):
    """log(exp(a) - exp(b)) with a >= b."""
    if b == NEG:
        return a
    m = max(a, b)
    return m + math.log(math.exp(a - m) - math.exp(b - m))


def coloc(d1, d2, p1=1e-4, p2=1e-4, p12=1e-5):
    """Return (n_snp, [PP.H0..PP.H4]) for two {rsid:(beta,se)} dicts."""
    common = sorted(set(d1) & set(d2))
    if not common:
        return 0, None
    la = [labf(*d1[rs]) for rs in common]
    lb = [labf(*d2[rs]) for rs in common]
    l1 = logsumexp(la)
    l2 = logsumexp(lb)
    l12 = logsumexp([x + y for x, y in zip(la, lb)])
    l3 = log_diff(l1 + l2, l12)          # s1*s2 - s12  ==  H3 term
    logH = [0.0, math.log(p1) + l1, math.log(p2) + l2, math.log(p1) + math.log(p2) + l3,
            math.log(p12) + l12]
    m = max(logH)
    tot = m + math.log(sum(math.exp(x - m) for x in logH))
    return len(common), [math.exp(x - tot) for x in logH]


def main():
    d1 = load(sys.argv[1]); d2 = load(sys.argv[2])
    n, pp = coloc(d1, d2)
    if pp is None:
        print("no common SNP"); return
    print(f"nSNP={n}  PP.H0={pp[0]:.3f} H1={pp[1]:.3f} H2={pp[2]:.3f} "
          f"H3={pp[3]:.3f} H4(共享)={pp[4]:.3f}")


if __name__ == "__main__":
    main()
