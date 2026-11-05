#!/usr/bin/env bash
# Run analysis/analyze.py on REAL reads: a public Illumina MiSeq ARTIC v3 SARS-CoV-2 amplicon run
# (98 amplicons, 2 tiled pools, ~27.7k read pairs) from the nf-core/test-datasets repository.
# Needs internet (raw.githubusercontent.com) and the conda env. Not part of CI or check_all.sh.
#   bash scripts/real_data_check.sh [workdir]
set -euo pipefail
cd "$(dirname "$0")/.."
W="${1:-$(mktemp -d)}"; mkdir -p "$W"
B=https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon
G=genome/NC_045512.2
fetch() { [ -s "$W/$(basename "$1")" ] || curl -fsSL -o "$W/$(basename "$1")" "$B/$1"; }
for f in illumina/amplicon/sample1_R1.fastq.gz illumina/amplicon/sample1_R2.fastq.gz \
         $G/GCF_009858895.2_ASM985889v3_genomic.200409.fna.gz $G/amplicon/nCoV-2019.artic.V3.bed \
         $G/amplicon/nCoV-2019.artic.V3.primer.fasta; do fetch "$f"; done
( cd "$W" && sha256sum -c --quiet - <<'S'
b7469350e3167dcdeab6a498984030b09af4bbf17e5b5a9e8f95dc3ee352b031  sample1_R1.fastq.gz
1fcbddf6dbb6d6508755477859161fa63be14246a85d4a4f75683da3f461e153  sample1_R2.fastq.gz
4a1c06f35d84707b60300d1250ef5f764dce80da90e85fe55a49a0e604ec8158  nCoV-2019.artic.V3.bed
4090454a089d8e7a2bbf01abc2dae5aad5d44ec82892d7d82ed136608f183a75  nCoV-2019.artic.V3.primer.fasta
d359795c5b11e9fc9309d5b3edee2dde1040781cd6a13b09aa7f4aaad8f07816  GCF_009858895.2_ASM985889v3_genomic.200409.fna.gz
S
)
gunzip -c "$W"/GCF_*.fna.gz > "$W/ref.fa"
# ARTIC primer scheme -> the inputs analyze.py expects: amplicon BED (outer primer edges) + primers named *_F / *_R
python - "$W" <<'P'
import collections, re, sys
w = sys.argv[1]
left, right = collections.defaultdict(list), collections.defaultdict(list)
for line in open(f"{w}/nCoV-2019.artic.V3.bed"):
    c, s, e, n = line.split()[:4]
    i, side = re.match(r"nCoV-2019_(\d+)_(LEFT|RIGHT)", n).groups()
    (left if side == "LEFT" else right)[int(i)].append((int(s), int(e)))
with open(f"{w}/amplicons.bed", "w") as bed:
    for i in sorted(left):
        bed.write(f"NC_045512.2\t{min(a for a, _ in left[i])}\t{max(b for _, b in right[i])}\tamp{i}\t0\t+\n")
seqs, name = {}, None
for line in open(f"{w}/nCoV-2019.artic.V3.primer.fasta"):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:]; seqs[name] = ""
    else:
        seqs[name] += line
with open(f"{w}/primers.fasta", "w") as fa:
    for n, s in seqs.items():
        i, alt, side = re.match(r"nCoV-2019_(\d+)(_alt\d*)?_(LEFT|RIGHT)", n).groups()
        fa.write(f">amp{i}{alt or ''}_{'F' if side == 'LEFT' else 'R'}\n{s}\n")
P
# this library starts reads at either end of an amplicon, hence --either-orientation
python analysis/analyze.py --reads "$W/sample1_R{1,2}.fastq.gz" --bed "$W/amplicons.bed" --primers "$W/primers.fasta" \
    --genome "$W/ref.fa" --outdir "$W/out" --either-orientation
python - "$W/out/sample1/sample1_summary.json" "$W/out/sample1/target_coverage.tsv" <<'P'
import csv, json, sys
s = json.load(open(sys.argv[1]))
cov = [float(r["mean_coverage"]) for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")]
assert s["n_targets"] == 98 == len(cov), s["n_targets"]
assert s["pct_retained"] > 70 and s["pct_reads_mapped"] > 99, (s["pct_retained"], s["pct_reads_mapped"])
assert min(cov) > 1, "an amplicon has no coverage"
print(f"ok: real data, {s['pct_retained']:.0f}% of pairs kept, fold-80 {s['fold_80']:.1f}, lowest amplicon {min(cov):.0f}x")
P
