"""Reproducible E. coli target sets from the tracked GFF (for design, simulation and tests).

    python scripts/make_ecoli_targets.py spread --n 100 --out design/targets/ecoli_spread_100.tsv
    python scripts/make_ecoli_targets.py tiled  --n 100 --out design/targets/ecoli_tiled_100.tsv

spread: N protein-coding genes (>= --min-length, no insertion-sequence or rhs repeats) evenly spaced by
        position along the genome, so panels are genome-wide, not clustered.
tiled:  N overlapping windows (--window, stepping --step) from --start, each with a 230-300 bp amplicon.
"""
import argparse
import re
import sys

CONTIG = "NC_000913.3"
REPEATS = re.compile(r"^(ins|rhs)", re.IGNORECASE)     # multi-copy families that cannot be specific


def read_genes(gff, min_length):
    """Unique protein-coding genes: [(name, contig, start, end)] sorted by position (1-based, inclusive)."""
    genes, seen = [], set()
    for line in open(gff):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9 or f[2] != "gene":
            continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        name = attrs.get("Name") or attrs.get("gene") or attrs.get("locus_tag")
        start, end = int(f[3]), int(f[4])
        if (not name or attrs.get("gene_biotype") != "protein_coding" or end - start + 1 < min_length
                or REPEATS.match(name) or name in seen):
            continue
        seen.add(name)
        genes.append((name, f[0], start, end))
    return sorted(genes, key=lambda g: (g[1], g[2]))


def spread(genes, n):
    if n > len(genes):
        raise ValueError(f"only {len(genes)} eligible genes, asked for {n}")
    return [genes[int(i * len(genes) / n)] for i in range(n)]


def tiled(n, start, window, step, amp_min=230, amp_max=300):
    return [(f"T{i + 1:03d}", CONTIG, start + i * step, start + i * step + window - 1, amp_min, amp_max) for i in range(n)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=["spread", "tiled"])
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--gff", default="common_refs/genomic.gff")
    ap.add_argument("--min-length", type=int, default=600)
    ap.add_argument("--start", type=int, default=100001)
    ap.add_argument("--window", type=int, default=300)
    ap.add_argument("--step", type=int, default=200)
    a = ap.parse_args()
    if a.mode == "spread":
        rows = [(n, c, s, e) for n, c, s, e in spread(read_genes(a.gff, a.min_length), a.n)]
        header = "name\treference\tstart\tend"
    else:
        rows = tiled(a.n, a.start, a.window, a.step)
        header = "name\treference\tstart\tend\tamplicon_min\tamplicon_max"
    with open(a.out, "w") as f:
        f.write(header + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))
    print(f"wrote {len(rows)} {a.mode} targets to {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
