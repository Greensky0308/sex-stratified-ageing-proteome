#!/usr/bin/env python3
"""Reverse-direction MR: lean mass -> SULT2A1 (test for reverse causation).

Instruments: distance-clumped genome-wide-significant ALM leads (Pei 2020),
excluding the SULT2A1 locus. Outcome: SULT2A1 aptamer (deCODE SomaScan,
genome-wide, combined sexes).
"""
import gzip, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOMA = os.path.join(ROOT, "data", "decode",
                    "Proteomics_SMP_PC0_9829_91_SULT2A1_SULT_2A1_10032022.txt.gz")
ALM = {"F": "GCST90000027.h.tsv.gz", "M": "GCST90000026.h.tsv.gz"}
SULTLOC = ("19", 46_000_000, 50_000_000)     # exclude cis-SULT2A1
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def alm_leads(acc, win=1_000_000, p_thr=5e-8):
    """Distance-clumped genome-wide-significant leads (one per window)."""
    hits = []
    with gzip.open(os.path.join(ROOT, "data", "sumstats", acc), "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            try:
                if float(c[ix["p_value"]]) >= p_thr:
                    continue
            except ValueError:
                continue
            ch = c[ix["chromosome"]].replace("chr", "")
            pos = int(c[ix["base_pair_location"]])
            hits.append((ch, pos, c[ix["effect_allele"]], c[ix["other_allele"]],
                         float(c[ix["beta"]]), float(c[ix["standard_error"]])))
    hits.sort(key=lambda h: (h[0], h[1]))
    leads, seen = [], {}
    for h in hits:
        key = (h[0], h[1] // win)
        if key not in seen:
            seen[key] = True
            leads.append(h)
    return leads


def scan_soma(targets):
    d = {}
    with gzip.open(SOMA, "rt") as fh:
        fh.readline()
        for line in fh:
            c = line.split("\t")
            if len(c) < 11:
                continue
            try:
                key = (c[0].replace("chr", ""), int(c[1]))
            except ValueError:
                continue
            if key in targets:
                try:
                    d[key] = (c[4], c[5], float(c[6]), float(c[9]))   # ea,oa,beta,se
                except ValueError:
                    continue
    return d


def harm(beta, ea_o, oa_o, ea_s, oa_s):
    if ea_s == ea_o and oa_s == oa_o:
        return beta
    if ea_s == oa_o and oa_s == ea_o:
        return -beta
    if ea_s == COMP.get(ea_o) and oa_s == COMP.get(oa_o):
        return beta
    if ea_s == COMP.get(oa_o) and oa_s == COMP.get(ea_o):
        return -beta
    return None


def main():
    per_sex = {}
    allt = {}
    for sex in ("F", "M"):
        leads = alm_leads(ALM[sex])
        leads = [h for h in leads
                 if not (h[0] == SULTLOC[0] and SULTLOC[1] <= h[1] <= SULTLOC[2])]
        per_sex[sex] = leads
        print(f"ALM {sex}: clumped leads (excl. SULT2A1 locus) = {len(leads)}")
        for h in leads:
            allt[(h[0], h[1])] = True
    soma = scan_soma(set(allt))
    print(f"SomaScan SNPs matched at ALM-lead positions = {len(soma)}")

    for sex in ("F", "M"):
        num = den = 0.0
        used = 0
        for (ch, pos, ea, oa, b, se) in per_sex[sex]:
            s = soma.get((ch, pos))
            if not s:
                continue
            s_ea, s_oa, s_b, s_se = s
            # exposure (ALM) per SomaScan allele -> Wald = b_SULT2A1 / b_ALM
            b_alm = harm(b, ea, oa, s_ea, s_oa)   # ALM beta per SomaScan EA
            if b_alm is None or abs(b_alm) < 1e-9:
                continue
            wz = s_b / b_alm
            wse = s_se / abs(b_alm)
            num += wz / wse ** 2
            den += 1 / wse ** 2
            used += 1
        if used:
            ivw = num / den
            se_ivw = math.sqrt(1 / den)
            print(f"  ALM->SULT2A1 {sex}: nSNP={used}  IVW beta={ivw:+.4f} "
                  f"se={se_ivw:.4f}  z={ivw/se_ivw:+.2f}")


if __name__ == "__main__":
    main()
