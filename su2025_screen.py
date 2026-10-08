#!/usr/bin/env python3
"""Full cis-pQTL sex-stratified MR screen in Su 2025 (TOPMed+LOS) A-LM / WB-LM.

Cross-cohort check: does SULT2A1 still rank top when the discovery screen is
re-run in an independent, non-UKB, sex-stratified cohort (even if underpowered)?
"""
import gzip, math, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KP4 = os.path.join(ROOT, "data", "kp4cd")
INST = os.path.join(ROOT, "data", "instruments_ukbppp_grch38.tsv")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def load_instruments():
    d = {}
    with open(INST) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            prot, rs, ch, pos, ea, bx = c[0], c[1], c[2], int(c[3]), c[4], float(c[5])
            if abs(bx) < 1e-9:
                continue
            d[(ch, pos)] = (prot, ea, bx)
    return d


def scan_su(fn, targets):
    """Return {pos: (a1,a2,eff,se,p,n)} for chr19 target positions."""
    d = {}
    path = os.path.join(KP4, fn)
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        i_mk, i_a1, i_a2 = ix["MarkerName"], ix["Allele1"], ix["Allele2"]
        i_ef, i_se, i_p = ix["Effect"], ix["StdErr"], ix["P-value"]
        try:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) <= i_p:
                    continue
                mk = c[i_mk]
                c1 = mk.find(":")
                c2 = mk.find(":", c1 + 1)
                if c1 < 0 or c2 < 0:
                    continue
                try:
                    key = (mk[:c1], int(mk[c1 + 1:c2]))
                except ValueError:
                    continue
                if key not in targets:
                    continue
                try:
                    d[key] = (c[i_a1].upper(), c[i_a2].upper(),
                              float(c[i_ef]), float(c[i_se]), float(c[i_p]))
                except ValueError:
                    continue
        except (EOFError, gzip.BadGzipFile):
            pass
    return d


def align(bx, ea_inst, a1, a2):
    """Instrument beta expressed per Su2025 allele a1."""
    if a1 == ea_inst:
        return bx
    if a2 == ea_inst:
        return -bx
    if a1 == COMP.get(ea_inst):
        return bx
    if a2 == COMP.get(ea_inst):
        return -bx
    return None


def main():
    inst = load_instruments()
    targets = set(inst.keys())
    print(f"instruments={len(inst)}")

    per = {"F": {}, "M": {}}
    for sex, fn in (("F", "su2025_applean_F.tsv.gz"), ("M", "su2025_applean_M.tsv.gz")):
        su = scan_su(fn, targets)
        print(f"  {fn}: matched target SNPs = {len(su)}", flush=True)
        for key, (a1, a2, eff, se, p) in su.items():
            if se <= 0:
                continue
            prot, ea, bx = inst[key]
            b = align(bx, ea, a1, a2)
            if b is None:
                continue
            per[sex][prot] = (eff / b, se / abs(b))

    res = []
    for prot in per["F"]:
        if prot not in per["M"]:
            continue
        bF, sF = per["F"][prot]
        bM, sM = per["M"][prot]
        z = (bF - bM) / math.sqrt(sF ** 2 + sM ** 2)
        res.append((prot, bF, sF, bM, sM, z))
    res.sort(key=lambda r: abs(r[5]), reverse=True)
    out = os.path.join(ROOT, "out", "su2025_screen_sexdiff.tsv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as fh:
        fh.write("protein\tbF\tseF\tbM\tseM\tzdiff\n")
        for r in res:
            fh.write(f"{r[0]}\t{r[1]:.4g}\t{r[2]:.4g}\t{r[3]:.4g}\t{r[4]:.4g}\t{r[5]:.3f}\n")
    print(f"paired proteins={len(res)} -> {out}")
    rank = next((i for i, r in enumerate(res, 1) if r[0] == "SULT2A1"), None)
    print("top 12 by |z|:")
    for i, r in enumerate(res[:12], 1):
        print(f"  {i:2d} {r[0]:12s} bF={r[1]:+.4f} bM={r[3]:+.4f} z={r[5]:+.2f}")
    if rank:
        r = res[rank - 1]
        print(f"SULT2A1 rank = {rank}/{len(res)}  z={r[5]:+.2f}  bF={r[1]:+.4f} bM={r[3]:+.4f}")
    else:
        print("SULT2A1 not in Su2025 paired set")


if __name__ == "__main__":
    main()
