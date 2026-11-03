"""Percentage of genome positions whose k-mer occurs more than once, for a range of k.

Vectorised with numpy (2-bit encoding, k <= 31). Windows containing non-ACGT bases
are skipped. With --canonical, a k-mer and its reverse complement count as the same.
"""
import argparse
import os

import numpy as np
from Bio import SeqIO

LUT = np.full(256, 255, dtype=np.uint8)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i
    LUT[b + 32] = i  # lowercase


def encode_genome(fasta_path):
    """Return a list of uint8 arrays (one per record), 0-3 for ACGT, 255 otherwise."""
    return [LUT[np.frombuffer(str(r.seq).encode(), dtype=np.uint8)]
            for r in SeqIO.parse(fasta_path, "fasta")]


def kmer_codes(arr, k, canonical):
    """2-bit packed k-mer codes for every valid window of one record."""
    n = len(arr) - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.uint64)
    a = arr.astype(np.uint64)
    fwd = np.zeros(n, dtype=np.uint64)
    rc = np.zeros(n, dtype=np.uint64)
    bad = np.zeros(n, dtype=bool)
    for j in range(k):
        col = a[j:j + n]
        bad |= col > 3
        c = np.where(col > 3, 0, col)
        fwd = (fwd << np.uint64(2)) | c
        if canonical:
            rc |= (np.uint64(3) - c) << np.uint64(2 * j)
    codes = np.minimum(fwd, rc) if canonical else fwd
    return codes[~bad]


def measure_uniqueness(fasta_path, start_k=12, end_k=25, canonical=False):
    records = encode_genome(fasta_path)
    print(f"Analyzing genome: {fasta_path} ({sum(len(r) for r in records)} bp)")
    print(f"{'K-mer':<8} | {'Total':<12} | {'Non-unique':<12} | {'% non-unique':<12}")
    ks, pcts = [], []
    for k in range(start_k, end_k + 1):
        codes = np.concatenate([kmer_codes(r, k, canonical) for r in records])
        total = len(codes)
        if total == 0:
            continue
        _, counts = np.unique(codes, return_counts=True)
        non_unique = int(counts[counts > 1].sum())  # positions whose k-mer repeats
        pct = 100.0 * non_unique / total
        ks.append(k)
        pcts.append(pct)
        print(f"{k:<8} | {total:<12} | {non_unique:<12} | {pct:.2f}%")
    return ks, pcts


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("fasta")
    ap.add_argument("--start-k", type=int, default=12)
    ap.add_argument("--end-k", type=int, default=25)
    ap.add_argument("--canonical", action="store_true")
    ap.add_argument("--plot", default="kmer_uniqueness_plot.png")
    args = ap.parse_args()
    ks, pcts = measure_uniqueness(args.fasta, args.start_k, args.end_k, args.canonical)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.plot(ks, pcts, marker="o")
    plt.title(f"Genomic non-uniqueness vs k-mer length\n{os.path.basename(args.fasta)}")
    plt.xlabel("k-mer length (bp)")
    plt.ylabel("% positions with a repeated k-mer")
    plt.grid(True)
    plt.xticks(ks)
    plt.ylim(0, 100)
    plt.savefig(args.plot)
    print(f"Plot saved to {args.plot}")


if __name__ == "__main__":
    main()
