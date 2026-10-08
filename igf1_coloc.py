#!/usr/bin/env python3
"""Colocalisation: SULT2A1 cis-pQTL (SomaScan, GRCh38) vs IGF-1 (Neale 30770).

Tests whether the SULT2A1 sex-opposite effect on IGF-1 reflects a shared causal
variant at the locus or just LD. IGF-1 regional stats are pulled as a BGZF tail.
"""
import gzip, io, os, subprocess, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coloc_abf import coloc                                   # noqa: E402
from sult2a1_crossplatform import read_somascan, SOMA         # noqa: E402
from fetch_neale_snp import size_of, PROXY, BGZF, CACHE, CURL       # noqa: E402

from pyliftover import LiftOver                               # noqa: E402

B = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs"
LO37, HI37 = 47_600_000, 48_900_000      # GRCh37 window around SULT2A1
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def fetch_region(url, lo, hi, frac=0.14, tries=4):
    S = size_of(url)
    start = int(S * (1 - frac))
    os.makedirs(CACHE, exist_ok=True)
    tmp = os.path.join(CACHE, "_igf1_tail.bgz")
    for _ in range(tries):
        subprocess.run(CURL + ["-r", f"{start}-{S-1}", url, "-o", tmp], capture_output=True)
        if not os.path.exists(tmp):
            continue
        data = open(tmp, "rb").read()
        i = data.find(BGZF)
        if i < 0:
            continue
        d = {}
        try:
            with gzip.open(io.BytesIO(data[i:]), "rt") as fh:
                for line in fh:
                    if not line.startswith("19:"):
                        continue
                    c = line.split("\t")
                    mk = c[0].split(":")
                    p = int(mk[1])
                    if lo <= p <= hi:
                        d[p] = (mk[2], mk[3], float(c[7]), float(c[8]))
        except (EOFError, gzip.BadGzipFile, UnicodeDecodeError):
            pass
        if d:
            return d
    return {}


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
    ss = read_somascan(SOMA)                     # {rsid: pos(GRCh38), ea, oa, beta, se}
    bypos = {}
    for rs, r in ss.items():
        bypos.setdefault(r["pos"], []).append(r)
    lo = LiftOver("hg19", "hg38")

    for sex, field in (("F", "30770"), ("M", "30770")):
        url = f"{B}/{field}_irnt.gwas.imputed_v3.{'female' if sex=='F' else 'male'}.varorder.tsv.bgz"
        reg = fetch_region(url, LO37, HI37)
        print(f"\nIGF-1 {sex}: region SNPs (GRCh37) = {len(reg)}")
        table = {}
        for p37, (ref, alt, b, se) in reg.items():
            res = lo.convert_coordinate("chr19", p37)
            if not res:
                continue
            p38 = res[0][1]
            if p38 not in bypos:
                continue
            for r in bypos[p38]:
                bh = harm(b, ref, alt, r["ea"], r["oa"])
                if bh is None:
                    continue
                table[r["pos"]] = (bh, se)
                break
        if not table:
            print(f"  no overlapping SNPs -> skip")
            continue
        # window-stability check: coloc across shrinking windows around the gene
        # (guard against a neighbouring independent signal inflating H4)
        MID38 = 47_871_693
        for w in (1_300_000, 500_000, 250_000, 100_000, 50_000):
            sub = {pos: (b, se) for pos, (b, se) in table.items()
                   if abs(pos - MID38) <= w}
            ss_sub = {pos: (bypos[pos][0]["beta"], bypos[pos][0]["se"])
                      for pos in sub}
            if len(sub) < 50:
                print(f"  +/-{w//1000}kb: too few SNPs ({len(sub)})")
                continue
            n, pp = coloc(ss_sub, sub)
            print(f"  +/-{w//1000:>4d}kb  n={n:5d}  H3={pp[3]:.3f}  H4(shared)={pp[4]:.3f}")


if __name__ == "__main__":
    main()
