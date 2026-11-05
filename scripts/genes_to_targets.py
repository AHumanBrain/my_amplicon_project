"""Convert a list of gene names + a GFF3 into a targets TSV (one whole-gene region per row).

    python scripts/genes_to_targets.py --gff genomic.gff --genes genes.txt > targets.tsv

Matches Name=, gene= or locus_tag= attributes of 'gene' features. Positions are 1-based inclusive.
"""
import argparse
import re
import sys


def parse_gff(path):
    attr = re.compile(r"(?:Name|gene|locus_tag)=([^;]+)")
    genes = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'gene':
                continue
            for name in attr.findall(p[8]):
                genes.setdefault(name, (p[0], int(p[3]), int(p[4])))
    return genes


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--gff', required=True)
    ap.add_argument('--genes', required=True, help="Text file, one gene per line ('#' comments allowed)")
    a = ap.parse_args()
    genes = parse_gff(a.gff)
    print("name\treference\tstart\tend")
    missing = []
    for line in open(a.genes):
        g = line.strip()
        if not g or g.startswith('#'):
            continue
        if g in genes:
            print(f"{g}\t{genes[g][0]}\t{genes[g][1]}\t{genes[g][2]}")
        else:
            missing.append(g)
    if missing:
        print(f"warning: not found in GFF: {', '.join(missing)}", file=sys.stderr)


if __name__ == "__main__":
    main()
