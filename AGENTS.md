# AGENTS.md: guide for AI assistants working in this repo

A platform for **custom-panel multiplex PCR**: design primers for user-defined DNA or RNA targets, pack them into pools,
order them, analyse the sequencing, and feed coverage back into pool rebalancing. The user decides *what* to target; this
tool designs, screens and packs. Read `README.md` for detail; this file is the operating manual.

## Guardrails (follow these without being asked)
1. **Never order, pay, or contact a vendor.** Produce `<prefix>_order.csv` and stop; a human approves and orders.
2. **Show the design before it is ordered.** Run the viewer (below) and report pools, failed targets, flags and
   tiling gaps. Do not call a design "ready" if `results.failed` or `results.tiling_gaps` in the manifest is non-empty without saying so.
3. **Everything here is in silico.** Nothing is wet-lab validated. Say "predicted", not "works". The analysis was validated on simulated reads plus one real public dataset (ARTIC v3 SARS-CoV-2, `scripts/real_data_check.sh`); nothing human.
4. **Do not silently change thresholds, seeds or defaults.** If you change an option, say which and why. Same inputs + same `--seed` must reproduce identical output; keep it that way.
5. **No push to `main`, no merge, no force-push, no deleting data** unless the human explicitly asks in this session.
6. **Report failures plainly** (failed targets, unsatisfiable constraints, tests that fail). Never skip or weaken a test to get green.

## Workflow
```
 targets.tsv + reference FASTA(s)
        │  python design.py --project panel.toml          (or flags; CLI overrides the project file)
        ▼
 <prefix>_order.csv, _pool_N.{csv,bed}, _pool_N_trimming.fasta, _manifest.json, [_tiling_gaps.tsv, _failed_targets.txt]
        │  python scripts/design_viewer.py --prefix <prefix> --reference <fasta> --outdir <dir>
        ▼                                                  HUMAN REVIEW + APPROVAL, then the human orders
 wet lab (RT if RNA) → pooled PCR → Illumina sequencing → FASTQ
        │  python analysis/analyze.py --reads "x_R{1,2}.fastq.gz" --bed <prefix>_pool_N.bed \
        │      --primers <prefix>_pool_N_trimming.fasta --genome <reference.fa> --outdir results   (one run per pool)
        ▼
 results/<sample>/: sorted BAM, target_coverage.tsv, <sample>_balancing_feedback.json, _target_summary.tsv, _uniformity.png, _summary.json
        │  python design.py --feedback <json> --targets targets.tsv --project panel.toml
        ▼
 <prefix>_rebalancing_recipe.csv (+ redesign of targets flagged REDESIGN)
```
Start a new project by copying `panel.example.toml`. Targets: see the table in README section 2.

## Targets table (TSV, header required)
Columns: `name` (required, unique) · `reference start end` (1-based inclusive region on a contig/transcript in `--reference`) **or** `sequence` ·
`strand` (`+`/`-`: sense strand; design runs on the sense strand, forward primer = sense, reverse = antisense = RT primer for stranded RNA) ·
`junctions` (first base AFTER each splice join; a primer must overlap one) · `span` (amplicon must contain) · `avoid` (no primer overlaps) ·
`amplicon_min/max` · `max_hits`. Positions are 1-based within the target as given. A primer cannot overlap a `span` range.
For RNA: design on the **transcript** and pass `--specificity-fasta transcripts.fa genome.fa` so genomic amplification shows in `predicted_products`.

## Layouts
`--layout auto` picks **tiled** if target windows or best amplicons overlap, else **sparse**.
- **tiled**: at most 2 pools by position (A,C,.. | B,D,..); the Monte Carlo search only chooses primer pairs inside each pool.
  Cross-target products are blocked only below the max amplicon size (A and C share a pool ~2 amplicons apart).
  Check `tiling_gaps`: neighbours overlapping by < 2 primer lengths leave bases uncovered after primer trimming.
- **sparse**: fewest compatible pools (Monte Carlo first-fit); facing primers of different targets closer than 1000 bp cannot share a pool.
- **chain** (`--layout chain`, never chosen by `auto`): each target row is a REGION to cover completely (reference/start/end only). Sliding windows
  (`--chain-step`) supply candidate amplicons; `core/chain.py` chains them left to right with a guaranteed insert overlap (`--chain-overlap`),
  alternating 2 pools, backtracking when a choice dead-ends. Uncovered stretches are listed in `failed`. SARS-CoV-2 whole genome: 99.5% covered
  (ARTIC 99.6%), 114 amplicons in 2 pools, but a full run takes about 5 minutes (dimer thermodynamics against growing pools).

## Reading a result
`<prefix>_manifest.json` has inputs (sha256), options, layout, per-pool targets, `failed`, `tiling_gaps`, outputs and `next_steps`.
Per-pair CSV columns of note: `specificity_hits` (exact hits per database), `predicted_products` (in-silico PCR from exact hits).
Viewer near-matches (≥16 aligned bp, ≤1 mismatch, exact 3' end) are informational and not enforced by the designer.

## Verify before you say "done"
```bash
python -m pytest test_design.py -q                   # fast; BLAST is mocked; the end-to-end analysis test skips without bwa/samtools
python -m pyflakes design.py core analysis scripts
bash scripts/check_all.sh                            # everything with the REAL tools (conda env), ~1 min; run before any merge
```
Automatic CI (`.github/workflows/ci.yml`) is only lint + unit tests. `real-tools.yml` runs `check_all.sh` on demand. Change design logic ⇒ add a test;
change the analysis ⇒ run `check_all.sh`, a unit test alone is not enough. `research/` is not linted or tested.

## Map

`design.py` (design, packing, feedback, project/manifest) · `core/biophysics.py` (primer3 thermodynamics, KAPA conditions) · `core/uniqueness.py` (genome k-mer copy counts for the repeat screen) ·
`analysis/analyze.py` (+ `metrics.py`, `balance.py`, `plot_uniformity.py`, `generate_simulated_fastq.py`) ·
`scripts/` (`design_viewer.py`, `check_all.sh`, `real_data_check.sh` (real reads, needs internet), `genes_to_targets.py`, `make_*_example/fixture.py`) · `examples/sars2/` (SARS-CoV-2 reference + ARTIC v3 scheme) · `research/` (`plexity_study/`, `igv_showcase/`, `benchmark_sars2/`: unmaintained) · `panel.example.toml` · `test_design.py`.

## Decisions already made (do not reopen without the human)
SARS-CoV-2 (30 kb, `examples/sars2/`, `design/targets/sars2_*.tsv`) is the development and test genome: fast, relevant, and the real-read check uses it; E. coli is for scale tests only ·
primer size, GC and Tm are options (defaults 19-21 nt, 40-60%, Tm 59-61 are unchanged; they cannot design a 38% GC genome, use a `[primers]` section) · Tails are ON by default (truncated Illumina P5/P7; reads start at the gene-specific primer, so trimming uses specific primers only) ·
hairpin clamp deleted · no Docker · analysis is cutadapt + bwa + samtools + Python (Nextflow, GATK, Picard, FastQC, MultiQC removed; Picard parity verified) ·
no variant calling (amplicon reads are >99% start-duplicates: no duplicate-aware tools) · plexity study lives in `research/` · merge to `main` is the human's call.

## Open items (say so if relevant, do not paper over)
Tiling does not optimise neighbour overlap (it reports gaps) · strand specificity is decided at the RT step and cannot be verified from amplicon reads ·
dimer screening uses gene-specific sequence, not the tailed oligo · packing is a heuristic with no optimality guarantee · near-match search is informational ·
cross-product distance (1000 bp) is a chosen default, not a validated one · our thermodynamic metrics do not predict real amplicon coverage (benchmark_sars2: |rho| <= 0.19) · the repeat screen (`--uniqueness-screen`) is a conservative pre-filter (catches about half of repeat failures on human, no false rejections seen on 79 regions) and is off by default for genomes over ~500 Mb · targets must be chosen for specificity: random or repeat-derived regions are out of scope, the tool fails them rather than designing non-specific primers.
The in-silico plexity study (README section 4, code in `research/`) rests on stated assumptions: hot-start is modelled as dimers extending at 55 °C instead of 37 °C, "risky" = occupancy ≥ 0.42% (calibrated to the tool's own limit, no sensitivity run), the plain baseline is 250 nM per primer at any plexity, and the uniformity model is uncalibrated. Quote it as a model, never as a result.

## Ask the human when
the targets or reference are ambiguous · constraints are unsatisfiable · a design has failed targets or gaps · anything would leave the machine (ordering, publishing).
