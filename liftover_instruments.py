#!/usr/bin/env python3
"""LiftOver UKB-PPP cis instruments GRCh37 -> GRCh38 with pyliftover (local chain).

Input : data/instruments_ukbppp_grch37.tsv   (protein, rsid, chr37, pos37, ea, bx)
Output: data/instruments_ukbppp_grch38.tsv (protein, rsid, chr38, pos38, ea, bx)
"""
import os

from pyliftover import LiftOver

ROOT = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(ROOT, "data", "instruments_ukbppp_grch37.tsv")
DST = os.path.join(ROOT, "data", "instruments_ukbppp_grch38.tsv")


def main():
    lo = LiftOver("hg19", "hg38")
    n = bad = 0
    with open(SRC) as fh, open(DST, "w") as out:
        out.write("protein\trsid\tchr38\tpos38\tea\tbx\n")
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            prot, rs, ch37, pos37, ea, bx = c[0], c[1], c[2], int(c[3]), c[4], c[5]
            res = lo.convert_coordinate(f"chr{ch37}", pos37)
            if not res:
                bad += 1
                continue
            ch38 = res[0][0].replace("chr", "")
            pos38 = res[0][1]
            out.write(f"{prot}\t{rs}\t{ch38}\t{pos38}\t{ea}\t{bx}\n")
            n += 1
    print(f"lifted={n} failed={bad} -> {DST}")


if __name__ == "__main__":
    main()
