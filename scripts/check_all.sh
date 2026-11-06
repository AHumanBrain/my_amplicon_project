#!/usr/bin/env bash
# Full local check with the REAL tools (BLAST+, bwa, samtools, cutadapt). About 2 minutes.
#   micromamba activate amplicon_pipeline && bash scripts/check_all.sh
# Run it before merging. The automatic CI only does the fast part (lint + unit tests).
# Scratch files go to $CHECK_WORKDIR (default: a new temp dir that is left in place, nothing is deleted).
set -euo pipefail
cd "$(dirname "$0")/.."
for t in python blastn makeblastdb bwa samtools cutadapt; do
  command -v "$t" >/dev/null || { echo "missing: $t (activate the conda env from environment.yml)"; exit 1; }
done
export WORK="${CHECK_WORKDIR:-$(mktemp -d)}"
mkdir -p "$WORK"
echo "scratch dir: $WORK"
step() { echo; echo "### $*"; }

step "lint (undefined names are errors)"
python -m pyflakes design.py core analysis scripts

step "unit tests (includes the end-to-end analysis test: real cutadapt + bwa + samtools)"
python -m pytest test_design.py -q

# SARS-CoV-2 (30 kb) is the development genome: real BLAST in seconds. E. coli is used only for the scale check below.
SARS=examples/sars2/NC_045512.2.fa
# primer limits for an AT-rich genome (the tool's defaults cannot place primers across much of SARS-CoV-2)
RELAXED="--primer-min-size 20 --primer-max-size 30 --primer-min-gc 30 --primer-min-tm 57 --primer-max-tm 64"

step "sparse design on SARS-CoV-2: its ORFs, real BLAST, manifest"
python design.py --targets design/targets/sars2_genes.tsv --reference $SARS --blast-dir "$WORK/db" --output-prefix "$WORK/orfs" \
  --amplicon-min 100 --amplicon-max 250 $RELAXED 2>&1 | grep -E "Layout|targets ->|Failed" | sed 's/.*- //' || true
python - <<'PY'
import json, os
m = json.load(open(os.environ["WORK"] + "/orfs_manifest.json"))
r = m["results"]
assert m["layout"] == "sparse" and r["targets_designed"] >= 8 and m["options"]["primer_min_size"] == 20, r
print("ok:", r["targets_designed"], "of", r["targets_requested"], "ORFs designed;", "failed:", r["failed"] or "none")
PY

step "tiled design from the example project (first 30 of 99 tiles) -> two alternating pools, then the viewer"
head -31 design/targets/sars2_tiles.tsv > "$WORK/tiles30.tsv"
python design.py --project panel.example.toml --targets "$WORK/tiles30.tsv" --blast-dir "$WORK/db" --output-prefix "$WORK/t30" \
  --monte-carlo-iters 3 2>&1 | grep -E "Layout|targets ->|Tiling" | sed 's/.*- //' | cut -c1-160 || true
python - <<'PY'
import csv, glob, json, os
td = os.environ["WORK"] + "/t30"
m = json.load(open(td + "_manifest.json"))
assert m["layout"] == "tiled" and m["results"]["n_pools"] == 2 and m["results"]["targets_designed"] >= 25, m["results"]
pool = {r["target_id"]: int(r["pool"]) for f in glob.glob(td + "_pool_*.csv") for r in csv.DictReader(open(f))}
at = lambda i: pool.get(f"W{i:03d}")
assert all(at(i) != at(i + 1) for i in range(1, 30) if at(i) and at(i + 1)), "neighbouring tiles must be in different pools"
assert all(at(i) == at(i + 2) for i in range(1, 29) if at(i) and at(i + 2)), "tiles two apart must share a pool"
print("ok:", m["results"]["targets_designed"], "of 30 tiles designed; failed:", len(m["results"]["failed"]), "; gaps:", len(m["results"]["tiling_gaps"]))
PY
python scripts/design_viewer.py --prefix "$WORK/t30" --reference $SARS --blast-dir "$WORK/db" --outdir "$WORK/t30_view" | head -1
test -s "$WORK/t30_view/specificity.html"

step "RNA junction example: a junction-spanning primer must remove the genomic product"
python scripts/make_rna_example.py --outdir "$WORK/rna" >/dev/null
python design.py --targets "$WORK/rna/targets.tsv" --reference "$WORK/rna/cdna.fa" \
  --specificity-fasta "$WORK/rna/cdna.fa" "$WORK/rna/genome.fa" --blast-dir "$WORK/rna/db" \
  --output-prefix "$WORK/rna/out" --num-candidates 20 >/dev/null 2>&1
python - <<'PY'
import csv, glob, os
rows = {r["target_id"]: r for f in glob.glob(os.environ["WORK"] + "/rna/out_pool_*.csv") for r in csv.DictReader(open(f))}
assert set(rows) == {"TX1_plain", "TX1_junction"}, sorted(rows)   # same window: the tiled layout uses two pools
assert "genome:" in rows["TX1_plain"]["predicted_products"], "an unconstrained design should show a genomic product"
assert "genome:" not in rows["TX1_junction"]["predicted_products"], "a junction-spanning design must not amplify genomic DNA"
print("ok: junction constraint removes the genomic product:", rows["TX1_junction"]["predicted_products"])
PY

step "50 genome-wide E. coli genes: design, simulate reads, analyse (the whole loop at scale)"
python design.py --targets design/targets/ecoli_spread_50.tsv --blast-dir "$WORK/db" \
  --output-prefix "$WORK/s50" --monte-carlo-iters 5 2>&1 | grep -E "Layout|targets ->" | sed 's/.*INFO - //'
python analysis/generate_simulated_fastq.py --genome common_refs/ecoli_genome.fna --bed "$WORK/s50_pool_1.bed" \
  --output-dir "$WORK/s50_reads" --sample-name S50 --coverage 100 --seed 1 >/dev/null
python analysis/analyze.py --reads "$WORK/s50_reads/S50_R{1,2}.fastq.gz" --genome common_refs/ecoli_genome.fna \
  --bed "$WORK/s50_pool_1.bed" --primers "$WORK/s50_pool_1_trimming.fasta" --outdir "$WORK/s50_res"
python - <<'PY'
import json, os
w = os.environ["WORK"]
s = json.load(open(w + "/s50_res/S50/S50_summary.json"))
n_pool1 = sum(1 for _ in open(w + "/s50_pool_1.bed"))
assert s["n_targets"] == n_pool1 >= 25, (s["n_targets"], n_pool1)
assert s["pct_retained"] > 99 and s["pct_reads_mapped"] > 99 and s["pct_on_target_bases"] > 95, s
assert s["uniformity_0.2x"] > 90, s["uniformity_0.2x"]
print(f"ok: {n_pool1} amplicons, mean {s['mean_target_coverage']:.0f}x, {s['uniformity_0.2x']:.0f}% above 0.2x mean")
PY

step "analysis from the command line on the small fixture (trim, align, coverage, feedback)"
python scripts/make_test_fixture.py --outdir "$WORK/fx" >/dev/null
python analysis/analyze.py --reads "$WORK/fx/FixtureSample_R{1,2}.fastq.gz" --genome "$WORK/fx/genome.fna" \
  --bed "$WORK/fx/targets.bed" --primers "$WORK/fx/primers_trimming.fasta" --outdir "$WORK/res"
for f in FixtureSample.sorted.bam FixtureSample_balancing_feedback.json FixtureSample_target_summary.tsv \
         FixtureSample_uniformity.png FixtureSample_summary.json target_coverage.tsv cutadapt.log; do
  test -s "$WORK/res/FixtureSample/$f" || { echo "missing or empty: $f"; exit 1; }
done

echo; echo "ALL CHECKS PASSED"
