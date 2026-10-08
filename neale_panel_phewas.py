#!/usr/bin/env python3
"""Sex-stratified PheWAS of the SULT2A1 cis-pQTL across a Neale UKB panel.

Fetches only the tail BGZF block containing chr19:48374950 from each trait file.
Harmonises to the SULT2A1-lowering allele A and reports F/M effects + sex-diff z.
"""
import math, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fetch_neale_snp import fetch_tail, parse   # noqa: E402

B = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs"
EA = "A"
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

PANEL = [
    # lean-mass axis (replication across BIA body composition)
    ("23101", "Whole-body fat-free mass"), ("23102", "Whole-body water mass"),
    ("23125", "Arm fat-free mass (L)"), ("23113", "Leg fat-free mass (R)"),
    ("23105", "Basal metabolic rate"),
    # body size / adiposity (specificity controls)
    ("23100", "Whole-body fat mass"), ("23124", "Arm fat mass (L)"),
    ("23112", "Leg fat mass (R)"), ("23104", "BMI"), ("21002", "Weight"),
    ("50", "Height"), ("3148", "Heel BMD"),
    # strength
    ("46", "Grip strength (L)"), ("47", "Grip strength (R)"),
    # hormonal / metabolic / renal / cognitive
    ("30850", "Testosterone"), ("30830", "SHBG"), ("30800", "Oestradiol"),
    ("30770", "IGF-1"), ("30710", "C-reactive protein"), ("30750", "HbA1c"),
    ("30650", "AST"), ("30760", "HDL cholesterol"), ("30780", "LDL cholesterol"),
    ("30870", "Urate"), ("30720", "Cystatin C"), ("30700", "Creatinine"),
    ("30600", "Albumin"), ("30840", "Total bilirubin"), ("30740", "Glucose"),
    ("30690", "Total cholesterol"), ("4079", "Diastolic BP"), ("4080", "Systolic BP"),
    ("20023", "Reaction time"),
]


def harmon(beta, ref, alt):
    if alt == EA or alt == COMP.get(EA):
        return beta
    if ref == EA or ref == COMP.get(EA):
        return -beta
    return None


def get(field, sex):
    for suf in ("", ".varorder"):
        url = f"{B}/{field}_irnt.gwas.imputed_v3.{sex}{suf}.tsv.bgz"
        ln = fetch_tail(url)
        if ln:
            d = parse(ln)
            b = harmon(d["beta"], d["ref"], d["alt"])
            if b is not None:
                return (b, d["se"], d["p"], d["n"])
    return None


def main():
    out = []
    print(f"{'field':>6s} {'trait':26s}{'bF':>10s}{'bM':>10s}{'zdiff':>8s}{'p':>9s}")
    for field, lab in PANEL:
        rF = get(field, "female")
        rM = get(field, "male")
        if not rF or not rM:
            print(f"{field:>6s} {lab:26s}  (missing)")
            continue
        bF, sF, pF, nF = rF
        bM, sM, pM, nM = rM
        z = (bF - bM) / math.sqrt(sF ** 2 + sM ** 2)
        p = math.erfc(abs(z) / math.sqrt(2))
        print(f"{field:>6s} {lab:26s}{bF:>+10.4f}{bM:>+10.4f}{z:>+8.2f}{p:>9.2g}"
              f"   pF={pF:.1g} pM={pM:.1g} nF={nF} nM={nM}")
        out.append((field, lab, bF, sF, bM, sM, z, p, pF, pM, nF, nM))
    dst = os.path.join(ROOT, "out", "sult2a1_neale_phewas.tsv")
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "w") as fh:
        fh.write("field\ttrait\tbF\tseF\tbM\tseM\tzdiff\tp\tpF\tpM\tnF\tnM\n")
        for r in out:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print("wrote", dst)


if __name__ == "__main__":
    main()
