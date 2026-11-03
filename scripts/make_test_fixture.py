"""Build a tiny, realistic test fixture (genome slice, BED, primers, simulated reads).

Used by CI and local smoke tests so no tracked file is ever overwritten:
    python scripts/make_test_fixture.py --outdir "$RUNNER_TEMP/fixture"
"""
import argparse
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "analysis"))
sys.path.insert(0, ROOT)
import generate_simulated_fastq as sim  # noqa: E402
from core.adapters import Adapters  # noqa: E402

AMPLICONS = [(1000, 1180), (3000, 3190), (5200, 5400)]
PRIMER_LEN = 20


def build(genome_path, outdir, slice_len=8000, coverage=50, seed=1):
    os.makedirs(outdir, exist_ok=True)
    rec = next(SeqIO.parse(genome_path, "fasta"))
    rec.description = ""
    rec.seq = rec.seq[:slice_len]
    genome = os.path.join(outdir, "genome.fna")
    SeqIO.write(rec, genome, "fasta")

    seq = str(rec.seq).upper()
    bed = os.path.join(outdir, "targets.bed")
    fasta = os.path.join(outdir, "primers_trimming.fasta")
    with open(bed, "w") as b, open(fasta, "w") as f:
        for i, (s, e) in enumerate(AMPLICONS, 1):
            name = f"t{i}_amplicon_Strategy_0_0"
            b.write(f"{rec.id}\t{s}\t{e}\t{name}\t0\t+\n")
            fwd = seq[s:s + PRIMER_LEN]
            rev = str(Seq(seq[e - PRIMER_LEN:e]).reverse_complement())
            f.write(f">t{i}_F\n{Adapters.FWD_P5_TRUNCATED}{fwd}\n")
            f.write(f">t{i}_R\n{Adapters.REV_P7_TRUNCATED}{rev}\n")

    sim.generate_fastq_records(sim.parse_bed(bed), {rec.id: rec}, coverage, 100,
                               os.path.join(outdir, "FixtureSample"), seed=seed)
    return outdir


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", default="common_refs/ecoli_genome.cleaned.fna")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    print("Fixture written to", build(a.genome, a.outdir))
