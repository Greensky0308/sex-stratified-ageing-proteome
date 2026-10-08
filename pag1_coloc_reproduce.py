#!/usr/bin/env python3
"""Reproduce the PAG1 colocalization reported for the brain-ageing layer.

PAG1 pQTL (UKB-PPP Olink, GRCh38) vs the five strongest sex-stratified imaging
phenotypes for PAG1.  Both files are GRCh38, so SNPs are matched by position and
harmonised to the pQTL effect allele.
"""
import gzip, os, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
from coloc_abf import coloc                                  # noqa: E402

PAG1 = os.path.join(ROOT, "data", "decode",
                    "GBR_UKB_OLINK_OID20108_PAG1_Phosphoprotein_associated_with_"
                    "glycosphingolipid_enriched_microdomains_1_adjAgeSexBatPC_"
                    "InvNorm_22122022.txt.gz")
LEAD = 81_134_351          # GRCh38, lead cis-pQTL by |z| in the local file
PHENOS = ["GCST90430212", "GCST90431261", "GCST90431262", "GCST90431979", "GCST90431614"]
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def read_pag1(lo, hi):
    d = {}
    with gzip.open(PAG1, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        i = {k: h.index(k) for k in ("Chrom", "Pos", "effectAllele", "otherAllele", "Beta", "SE")}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if c[i["Chrom"]] != "chr8":
                continue
            p = int(c[i["Pos"]])
            if not (lo <= p <= hi):
                continue
            try:
                d[p] = (c[i["effectAllele"]], c[i["otherAllele"]],
                        float(c[i["Beta"]]), float(c[i["SE"]]))
            except ValueError:
                continue
    return d


def read_brain(acc, lo, hi):
    d = {}
    path = os.path.join(ROOT, "data", "sumstats", f"{acc}.h.tsv.gz")
    with gzip.open(path, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        i = {k: h.index(k) for k in
             ("chromosome", "base_pair_location", "effect_allele", "other_allele",
              "beta", "standard_error")}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if c[i["chromosome"]] not in ("8", "chr8"):
                continue
            p = int(c[i["base_pair_location"]])
            if p > hi:
                continue
            if p < lo:
                continue
            try:
                d[p] = (c[i["effect_allele"]], c[i["other_allele"]],
                        float(c[i["beta"]]), float(c[i["standard_error"]]))
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
    for w in (250_000, 500_000, 1_000_000):
        lo, hi = LEAD - w, LEAD + w
        pag = read_pag1(lo, hi)
        pag_d = {p: (v[2], v[3]) for p, v in pag.items()}
        print(f"\n=== window +/-{w//1000} kb (lead {LEAD}) — pQTL SNPs in window: {len(pag)}")
        for acc in PHENOS:
            br = read_brain(acc, lo, hi)
            shared = {}
            for p, (ea_b, oa_b, b_b, se_b) in br.items():
                if p not in pag:
                    continue
                bh = harm(b_b, ea_b, oa_b, pag[p][0], pag[p][1])
                if bh is None:
                    continue
                shared[p] = (bh, se_b)
            if len(shared) < 50:
                print(f"  {acc}: too few shared SNPs ({len(shared)})")
                continue
            prot = {p: pag_d[p] for p in shared}
            n, pp = coloc(prot, shared)
            print(f"  {acc}: n={n:5d}  H3={pp[3]:.3f}  H4(shared)={pp[4]:.3f}")


if __name__ == "__main__":
    main()
