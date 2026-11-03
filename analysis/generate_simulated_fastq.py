"""Simulate paired-end amplicon reads from a BED of amplicon intervals.

Deterministic for a given --seed. Reads are streamed to disk.
"""
import argparse
import gzip
import os
import random
import zlib

from Bio import SeqIO
from Bio.Seq import Seq

BASES = "ACGT"


def load_genome(genome_path):
    print(f"Loading genome from {genome_path}...")
    return SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))


def parse_bed(bed_path):
    intervals = []
    with open(bed_path) as f:
        for line in f:
            if line.strip() and not line.startswith(("#", "track", "browser")):
                parts = line.rstrip("\n").split("\t")
                intervals.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
    return intervals


def mutate_sequence(sequence, rng, mutation_rate=0.01):
    """Introduce random SNVs; returns (new_sequence, [(pos, ref, alt), ...])."""
    seq = list(sequence)
    mutations = []
    for i, ref in enumerate(seq):
        if rng.random() < mutation_rate:
            alt = rng.choice([b for b in BASES if b != ref])
            seq[i] = alt
            mutations.append((i, ref, alt))
    return "".join(seq), mutations


def add_errors(read, rng, error_rate):
    if error_rate <= 0:
        return read
    out = list(read)
    for i, ref in enumerate(out):
        if rng.random() < error_rate:
            out[i] = rng.choice([b for b in BASES if b != ref])
    return "".join(out)


def stable_offset(name):
    """Per-target offset in [-2.5, 2.5). crc32 is stable across runs (hash() is not)."""
    return (zlib.crc32(name.encode()) % 5) - 2.5


def pairs_for_amplicon(name, coverage, rng, non_uniform, thermal_bias):
    current = coverage
    if thermal_bias is not None:
        target_opt = 62.0 + stable_offset(name)
        scaling = max(0.05, 1.0 - abs(thermal_bias - target_opt) * 0.2)
        current = int(coverage * scaling)
    if non_uniform:
        return int(current * rng.uniform(0.1, 2.0))
    return current


def generate_fastq_records(intervals, genome_dict, coverage, read_len, output_prefix,
                           non_uniform=False, thermal_bias=None, seed=0,
                           mutation_rate=0.01, error_rate=0.001):
    rng = random.Random(seed)
    qual = "I" * read_len
    r1_filename = f"{output_prefix}_R1.fastq.gz"
    r2_filename = f"{output_prefix}_R2.fastq.gz"
    total = 0

    with gzip.open(r1_filename, "wt") as f1, gzip.open(r2_filename, "wt") as f2:
        for chrom, start, end, name in intervals:
            if chrom not in genome_dict:
                print(f"Warning: {chrom} not in genome. Skipping {name}.")
                continue
            ref_seq = str(genome_dict[chrom].seq[start:end]).upper()
            if not ref_seq:
                print(f"Warning: empty interval for {name}. Skipping.")
                continue
            amplicon, muts = mutate_sequence(ref_seq, rng, mutation_rate)
            n_pairs = pairs_for_amplicon(name, coverage, rng, non_uniform, thermal_bias)
            print(f"   -> {n_pairs} pairs for {name} ({len(muts)} SNVs)")

            for _ in range(n_pairs):
                r1 = add_errors(amplicon[:read_len].ljust(read_len, "N"), rng, error_rate)
                tail = amplicon[-read_len:].rjust(read_len, "N")
                r2 = add_errors(str(Seq(tail).reverse_complement()), rng, error_rate)
                header = f"SimRead:{total}:{name}"
                f1.write(f"@{header} 1:N:0:1\n{r1}\n+\n{qual}\n")
                f2.write(f"@{header} 2:N:0:1\n{r2}\n+\n{qual}\n")
                total += 1

    print(f"Wrote {total} read pairs to {r1_filename} and {r2_filename}")
    return total


def main():
    parser = argparse.ArgumentParser(description="Generate simulated FASTQ data from BED intervals.")
    parser.add_argument("--genome", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--sample-name", default="Simulated_Sample")
    parser.add_argument("--coverage", type=int, default=100, help="Read pairs per amplicon")
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--non-uniform", action="store_true")
    parser.add_argument("--thermal-bias", type=float,
                        help="Simulated annealing temperature (C); coverage falls off away from ~62C")
    parser.add_argument("--mutation-rate", type=float, default=0.01, help="Per-base SNV rate per amplicon")
    parser.add_argument("--error-rate", type=float, default=0.001, help="Per-base sequencing error rate")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    generate_fastq_records(parse_bed(args.bed), load_genome(args.genome), args.coverage,
                           args.read_length, os.path.join(args.output_dir, args.sample_name),
                           non_uniform=args.non_uniform, thermal_bias=args.thermal_bias,
                           seed=args.seed, mutation_rate=args.mutation_rate,
                           error_rate=args.error_rate)
    print("Simulation complete.")


if __name__ == "__main__":
    main()
