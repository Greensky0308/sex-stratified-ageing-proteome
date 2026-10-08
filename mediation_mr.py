#!/usr/bin/env python3
"""Two-step mediation: SULT2A1 -> IGF-1 -> lean mass (sex-stratified).

Step A (PheWAS, per SULT2A1-lowering allele A): SULT2A1 -> IGF-1.
Step B (MR): IGF-1 (genome-wide leads, Neale 30770, excl. the SULT2A1 locus) -> ALM (Pei 2020).
Indirect effect = beta_A(SULT2A1->IGF1) * beta(IGF1->ALM); proportion = indirect / beta_A(SULT2A1->ALM).
"""
import gzip, math, os

from pyliftover import LiftOver

ROOT = os.path.dirname(os.path.abspath(__file__))
LO = LiftOver("hg19", "hg38")
IGF = {"F": "data/igf1/igf1_female.tsv.bgz", "M": "data/igf1/igf1_male.tsv.bgz"}
ALM = {"F": "GCST90000027.h.tsv.gz", "M": "GCST90000026.h.tsv.gz"}
SULTLOC = ("19", 46_000_000, 50_000_000)
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

# Step A: SULT2A1 (per allele A) -> IGF-1  [from out/sult2a1_neale_phewas_all.tsv]
A_SULTI_IGF1 = {"F": 0.0338, "M": 0.0104}
# Step A total: SULT2A1 (per allele A) -> ALM  [converted from DXA Wald x bx]
A_SULTI_ALM = {"F": 0.01385, "M": -0.01207}


def igf_leads(sex, win=1_000_000, p_thr=5e-8):
    hits = []
    with gzip.open(os.path.join(ROOT, IGF[sex]), "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        try:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                try:
                    if float(c[ix["pval"]]) >= p_thr:
                        continue
                    fmt = c[ix["variant"]].split(":")
                    ch, pos37, ref, alt = fmt[0], int(fmt[1]), fmt[2], fmt[3]
                    b, se = float(c[ix["beta"]]), float(c[ix["se"]])
                    res = LO.convert_coordinate("chr" + ch, pos37)
                    if not res:
                        continue
                    ch, pos = res[0][0].replace("chr", ""), res[0][1]
                except (ValueError, KeyError, IndexError):
                    continue
                if ch == SULTLOC[0] and SULTLOC[1] <= pos <= SULTLOC[2]:
                    continue
                hits.append((ch, pos, ref, alt, b, se))
        except (EOFError, gzip.BadGzipFile):
            pass
    hits.sort(key=lambda h: (h[0], h[1]))
    leads, seen = [], set()
    for h in hits:
        key = (h[0], h[1] // win)
        if key not in seen:
            seen.add(key); leads.append(h)
    return leads


def alm_at(acc, targets):
    d = {}
    with gzip.open(os.path.join(ROOT, "data", "sumstats", acc), "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            key = (c[ix["chromosome"]].replace("chr", ""), int(c[ix["base_pair_location"]]))
            if key in targets:
                try:
                    d[key] = (c[ix["effect_allele"]], c[ix["other_allele"]],
                              float(c[ix["beta"]]), float(c[ix["standard_error"]]))
                except ValueError:
                    continue
    return d


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
    for sex in ("F", "M"):
        p = os.path.join(ROOT, IGF[sex])
        if not os.path.exists(p):
            print(f"[skip] {IGF[sex]}"); continue
        leads = igf_leads(sex)
        print(f"\nIGF-1 {sex}: clumped GW leads (excl. SULT2A1 locus) = {len(leads)}")
        targets = {(h[0], h[1]) for h in leads}
        alm = alm_at(ALM[sex], targets)
        print(f"  matched in ALM = {len(alm)}")
        num = den = 0.0; used = 0
        for (ch, pos, ref, alt, b_igf, se_igf) in leads:
            if (ch, pos) not in alm or abs(b_igf) < 1e-9:
                continue
            ea, oa, b_a, se_a = alm[(ch, pos)]
            # ALM beta per IGF-1 ALT allele
            b_al = harm(b_a, ea, oa, alt, ref) if alt in (ea, oa, COMP.get(ea), COMP.get(oa)) else None
            if b_al is None:
                continue
            w = b_al / b_igf
            ws = se_a / abs(b_igf)
            num += w / ws ** 2; den += 1 / ws ** 2; used += 1
        if not used:
            print("  no usable instruments"); continue
        b_igf_alm = num / den
        se_b = math.sqrt(1 / den)
        z = b_igf_alm / se_b
        print(f"  IGF-1 -> ALM ({sex}): nSNP={used} beta={b_igf_alm:+.4f} se={se_b:.4f} z={z:+.2f} p={math.erfc(abs(z)/math.sqrt(2)):.2g}")
        # mediation
        indirect = A_SULTI_IGF1[sex] * b_igf_alm
        total = A_SULTI_ALM[sex]
        if abs(total) > 1e-9:
            print(f"  indirect (SULT2A1->IGF1->ALM) = {indirect:+.5f}; total = {total:+.5f}; proportion = {indirect/total:.2%}")


if __name__ == "__main__":
    main()
