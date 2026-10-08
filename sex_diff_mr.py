#!/usr/bin/env python3
"""Sex-differential causal-proteome analysis for brain-aging IDP phenotypes.

For each (protein, phenotype) we run a two-sample Wald-ratio MR separately in
females and males (using cis-pQTL instruments), then test whether the causal
effect differs between sexes.

Aggregation across phenotypes uses Stouffer's method (sum(z)/sqrt(n)); proteins
are reported after Benjamini-Hochberg FDR.

Inputs
------
--mr        TSV from brain_extract:  protein  rsid  accession  b  se  z
--phenlist  JSON list of [tag, accession, trait]  (trait ends ' - female'/' - male')
--instruments TSV: protein rsid chr pos ea beta [p ...]
--gene-coords JSON: gene -> [chr, start, end]   (GRCh38)
--outdir    output directory

Outputs
-------
all.tsv  : protein, n_phen, z, p, q            (all instruments)
cis.tsv  : same, restricted to cis instruments (SNP within +/-1 Mb of gene)
trans.tsv: same, trans instruments
"""
import argparse, json, math, os, re
from collections import defaultdict

import numpy as np
from scipy.stats import norm


def benjamini_hochberg(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        prev = min(prev, p[idx] * n / (i + 1))
        q[idx] = prev
    return q


def base_trait(t):
    return re.sub(r"\s*-\s*(female|male)\s*$", "", t, flags=re.I).strip()


def sex_of(t):
    return "F" if t.lower().rstrip().endswith("female") else "M"


def load_instruments(path):
    """protein -> (gene, chr, pos). Gene parsed from 'NAME (GENE.seqid...)' or bare symbol."""
    inst = {}
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 4:
                continue
            m = re.search(r"\(([A-Za-z0-9\-]+)\.", c[0])
            gene = m.group(1) if m else c[0].split()[0]
            try:
                inst[c[0]] = (gene, str(c[2]), int(c[3]))
            except ValueError:
                continue
    return inst


def is_cis(protein, inst, coords, window=1_000_000):
    if protein not in inst:
        return None
    gene, chrom, pos = inst[protein]
    if gene not in coords:
        return None
    gc, gstart, gend = coords[gene]
    return str(chrom) == str(gc) and abs(pos - (gstart + gend) // 2) < window


def aggregate(mr_rows, phen_map, max_z=None):
    """Return list of (protein, n_phen, z) combining sex-difference z across phenotypes.

    max_z: if set, per-phenotype sex-difference z with |z| > max_z is dropped before
    aggregation (QC for implausible per-phenotype estimates). None = no filter.
    """
    acc_base = {acc: base_trait(t) for acc, t in phen_map.items()}
    acc_sex = {acc: sex_of(t) for acc, t in phen_map.items()}
    per_protein = defaultdict(dict)
    for row in mr_rows:
        protein, _rsid, acc, b, se, _z = row
        if acc not in acc_base:
            continue
        per_protein[(protein, acc_base[acc])].setdefault(acc_sex[acc], (b, se))
    dropped = 0
    total = 0
    zs_by_protein = defaultdict(list)
    for (protein, _base), bysex in per_protein.items():
        if "F" in bysex and "M" in bysex:
            (bF, seF), (bM, seM) = bysex["F"], bysex["M"]
            if seF <= 0 or seM <= 0:
                continue
            z = (bF - bM) / math.sqrt(seF**2 + seM**2)
            total += 1
            if not math.isfinite(z) or (max_z is not None and abs(z) > max_z):
                dropped += 1
                continue
            zs_by_protein[protein].append(z)
    if max_z is not None:
        print(f"QC(max_z={max_z}): dropped {dropped}/{total} protein-phenotype estimates "
              f"({100.0*dropped/max(total,1):.2f}%)")
    out = []
    for protein, zs in zs_by_protein.items():
        if len(zs) == 0:
            continue
        z = sum(zs) / math.sqrt(len(zs))
        out.append((protein, len(zs), z))
    return out


def write_result(rows, path):
    if not rows:
        open(path, "w").write("protein\tn_phen\tz\tp\tq\n")
        return 0
    p = [2 * (1 - norm.cdf(abs(z))) for _, _, z in rows]
    q = benjamini_hochberg(p)
    order = np.argsort(q)
    with open(path, "w") as fh:
        fh.write("protein\tn_phen\tz\tp\tq\n")
        for i in order:
            prot, n, z = rows[i]
            fh.write(f"{prot}\t{n}\t{z:.4f}\t{p[i]:.3g}\t{q[i]:.3g}\n")
    return int((np.array(q) < 0.05).sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mr", required=True)
    ap.add_argument("--phenlist", required=True)
    ap.add_argument("--instruments", required=True)
    ap.add_argument("--gene-coords", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--max-z", type=float, default=None,
                    help="QC: drop per-phenotype sex-difference |z| above this (default: no filter)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    phen = json.load(open(args.phenlist))
    phen_map = {p[1]: p[2] for p in phen}
    inst = load_instruments(args.instruments)
    coords = json.load(open(args.gene_coords))

    mr_rows = []
    for line in open(args.mr):
        c = line.rstrip("\n").split("\t")
        if len(c) != 6:
            continue
        try:
            mr_rows.append((c[0], c[1], c[2], float(c[3]), float(c[4]), float(c[5])))
        except ValueError:
            continue
    print(f"MR rows: {len(mr_rows)} | phenotypes in list: {len(phen_map)} | instruments: {len(inst)}")

    all_rows = aggregate(mr_rows, phen_map, args.max_z)
    cis_rows, trans_rows = [], []
    for protein, n, z in all_rows:
        flag = is_cis(protein, inst, coords)
        (cis_rows if flag is True else trans_rows if flag is False else []).append((protein, n, z))

    n_all = write_result(all_rows, os.path.join(args.outdir, "all.tsv"))
    n_cis = write_result(cis_rows, os.path.join(args.outdir, "cis.tsv"))
    n_tr = write_result(trans_rows, os.path.join(args.outdir, "trans.tsv"))
    print(f"proteins: all={len(all_rows)} cis={len(cis_rows)} trans={len(trans_rows)} unresolved={len(all_rows)-len(cis_rows)-len(trans_rows)}")
    print(f"q<0.05:  all={n_all}  cis={n_cis}  trans={n_tr}")
    for label, rows in [("CIS", cis_rows), ("TRANS", trans_rows)]:
        top = sorted(rows, key=lambda r: -abs(r[2]))[:6]
        print(f"\n{label} top:")
        for prot, n, z in top:
            p = 2 * (1 - norm.cdf(abs(z)))
            print(f"  {prot[:44]:46s} n={n:3d} z={z:+.2f} p={p:.2e}")


if __name__ == "__main__":
    main()
