#!/usr/bin/env bash
# Run analysis/analyze.py on REAL reads: a public Illumina MiSeq ARTIC v3 SARS-CoV-2 amplicon run
# (98 amplicons, 2 tiled pools, ~27.7k read pairs) from the nf-core/test-datasets repository.
# Only the reads (7 MB) are downloaded (raw.githubusercontent.com); the SARS-CoV-2 reference and ARTIC files are in examples/sars2.
# Needs internet and the conda env. Not part of CI or check_all.sh.
#   bash scripts/real_data_check.sh [workdir]
set -euo pipefail
cd "$(dirname "$0")/.."
W="${1:-$(mktemp -d)}"; mkdir -p "$W"
B=https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon
for r in 1 2; do [ -s "$W/sample1_R$r.fastq.gz" ] || curl -fsSL -o "$W/sample1_R$r.fastq.gz" "$B/illumina/amplicon/sample1_R$r.fastq.gz"; done
( cd "$W" && sha256sum -c --quiet - <<'S'
b7469350e3167dcdeab6a498984030b09af4bbf17e5b5a9e8f95dc3ee352b031  sample1_R1.fastq.gz
1fcbddf6dbb6d6508755477859161fa63be14246a85d4a4f75683da3f461e153  sample1_R2.fastq.gz
S
)
# reference, amplicon BED and primers are committed in examples/sars2 (see its README for how they were derived)
E=examples/sars2
# this library starts reads at either end of an amplicon, hence --either-orientation
python analysis/analyze.py --reads "$W/sample1_R{1,2}.fastq.gz" --bed $E/artic_v3_amplicons.bed --primers $E/artic_v3_primers.fasta \
    --genome $E/NC_045512.2.fa --outdir "$W/out" --either-orientation
python - "$W/out/sample1/sample1_summary.json" "$W/out/sample1/target_coverage.tsv" <<'P'
import csv, json, sys
s = json.load(open(sys.argv[1]))
cov = [float(r["mean_coverage"]) for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")]
assert s["n_targets"] == 98 == len(cov), s["n_targets"]
assert s["pct_retained"] > 70 and s["pct_reads_mapped"] > 99, (s["pct_retained"], s["pct_reads_mapped"])
assert min(cov) > 1, "an amplicon has no coverage"
print(f"ok: real data, {s['pct_retained']:.0f}% of pairs kept, fold-80 {s['fold_80']:.1f}, lowest amplicon {min(cov):.0f}x")
P
