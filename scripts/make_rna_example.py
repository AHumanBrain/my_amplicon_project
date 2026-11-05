"""Build a synthetic 3-exon gene, its spliced cDNA and a two-row targets table.

    python scripts/make_rna_example.py --outdir examples/rna_junction

Writes genome.fa (with 800 bp introns), cdna.fa (650 nt), targets.tsv. Use it to see how a
junction constraint removes genomic-DNA amplification:

    python design.py --targets examples/rna_junction/targets.tsv --reference examples/rna_junction/cdna.fa \\
        --specificity-fasta examples/rna_junction/cdna.fa examples/rna_junction/genome.fa
"""
import argparse
import os
import random

EXON_LENGTHS = (200, 200, 250)
INTRON_LENGTH = 800


def build(outdir, seed=11):
    rng = random.Random(seed)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    exons = [rand(n) for n in EXON_LENGTHS]
    introns = [rand(INTRON_LENGTH) for _ in EXON_LENGTHS[:-1]]
    genome = rand(1500) + exons[0] + introns[0] + exons[1] + introns[1] + exons[2] + rand(1500)
    cdna = "".join(exons)
    j1, j2 = EXON_LENGTHS[0] + 1, EXON_LENGTHS[0] + EXON_LENGTHS[1] + 1  # first base after each join
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "genome.fa"), "w") as f:
        f.write(f">chrG\n{genome}\n")
    with open(os.path.join(outdir, "cdna.fa"), "w") as f:
        f.write(f">TX1\n{cdna}\n")
    with open(os.path.join(outdir, "targets.tsv"), "w") as f:
        f.write("name\treference\tstart\tend\tjunctions\tamplicon_min\tamplicon_max\n"
                f"TX1_junction\tTX1\t1\t{len(cdna)}\t{j1},{j2}\t70\t150\n"
                f"TX1_plain\tTX1\t1\t{len(cdna)}\t\t70\t150\n")
    return outdir


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", required=True)
    print("Example written to", build(ap.parse_args().outdir))
