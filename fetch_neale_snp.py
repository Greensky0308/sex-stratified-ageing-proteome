#!/usr/bin/env python3
"""Fetch the SULT2A1 lead SNP from a remote BGZF Neale file using HTTP Range.

The Neale additive-tsvs are position-sorted BGZF; instead of downloading the
whole ~500 MB file we grab the tail (chr19 sits ~89% through the genome),
align to a BGZF block boundary, and decompress from there.
"""
import io, gzip, os, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, ".cache")
PROXY = os.environ.get("HTTPS_PROXY")
CURL = ["curl", "-s"] + (["-x", PROXY] if PROXY else [])
BGZF = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
TARGET = "19:48374950:"


def size_of(url):
    p = subprocess.run(CURL + ["-I", url],
                       capture_output=True, text=True)
    for ln in p.stdout.splitlines():
        if ln.lower().startswith("content-length:"):
            return int(ln.split(":", 1)[1])
    return None


def fetch_tail(url, frac=0.14, tries=4):
    """One ranged GET of the tail; tolerate a truncated transfer.

    chr19:48374950 sits ~89% through the genome, so a tail from ~86% will
    contain it even if the transfer is cut short.
    """
    S = size_of(url)
    if not S:
        return None
    start = int(S * (1 - frac))
    os.makedirs(CACHE, exist_ok=True)
    tmp = os.path.join(CACHE, "_tail.bgz")
    for _ in range(tries):
        subprocess.run(CURL + ["-r", f"{start}-{S-1}", url, "-o", tmp], capture_output=True)
        if not os.path.exists(tmp) or os.path.getsize(tmp) < 50_000:
            continue
        data = open(tmp, "rb").read()
        i = data.find(BGZF)
        if i < 0:
            continue
        try:
            with gzip.open(io.BytesIO(data[i:]), "rt") as fh:
                for line in fh:
                    if line.startswith(TARGET):
                        return line.rstrip("\n")
        except (EOFError, gzip.BadGzipFile, UnicodeDecodeError):
            pass
    return None


def parse(line):
    c = line.split("\t")
    ref, alt = c[0].split(":")[2], c[0].split(":")[3]
    return dict(ref=ref, alt=alt, beta=float(c[7]), se=float(c[8]),
                p=float(c[10]), n=int(c[4]))


if __name__ == "__main__":
    url = sys.argv[1]
    ln = fetch_tail(url)
    print(parse(ln) if ln else "NOT FOUND", "|", os.path.basename(url))
