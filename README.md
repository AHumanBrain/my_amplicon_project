# Multiplex Amplicon Toolkit

[![CI](https://github.com/AHumanBrain/my_amplicon_project/actions/workflows/ci.yml/badge.svg)](https://github.com/AHumanBrain/my_amplicon_project/actions/workflows/ci.yml)

Custom-panel multiplex PCR for DNA or RNA targets. You decide what to target; the tool designs the primers, screens them, packs them into pools, and analyses the sequencing so coverage feeds back into pool rebalancing.

1. **Design** primers (primer3 + BLAST+ + dimer screening) for regions or sequences you specify
2. **Pack** them into the fewest compatible pools (at most two for tiled designs)
3. **Analyse** amplicon sequencing and turn per-amplicon coverage into a rebalancing recipe

Primers are ordinary oligos with truncated Illumina tails (P5 forward, P7 reverse), orderable from any vendor. No custom chemistry. Everything here is **in silico**: nothing is wet-lab validated.

```mermaid
%%{init: {"theme":"base","themeVariables":{"fontFamily":"Inter, system-ui, -apple-system, Segoe UI, sans-serif","fontSize":"15px","background":"#fcfcfb","primaryColor":"#cde2fb","primaryBorderColor":"#2a78d6","primaryTextColor":"#0b0b0b","secondaryColor":"#d3f1e5","tertiaryColor":"#fde1d6","lineColor":"#52514e","textColor":"#0b0b0b","clusterBkg":"#fcfcfb","clusterBorder":"#e4e3df","edgeLabelBackground":"#fcfcfb"},"flowchart":{"curve":"basis","nodeSpacing":30,"rankSpacing":42,"padding":12}}}%%
flowchart LR
    T["Targets TSV<br/>+ reference"]:::step --> D["design.py<br/>primers · pools"]:::step
    D --> O["order.csv<br/>pool BEDs<br/>manifest"]:::out
    O -. "order · PCR · sequence" .-> R["FASTQ"]:::plain
    R --> A["analyze.py<br/>trim · align · coverage"]:::step
    O --> A
    A --> F["report<br/>feedback JSON"]:::out
    F -- "rebalance · redesign" --> D
    classDef step fill:#cde2fb,stroke:#2a78d6,stroke-width:1.6px,color:#0b0b0b
    classDef out fill:#d3f1e5,stroke:#1baf7a,stroke-width:1.6px,color:#0b0b0b
    classDef warn fill:#fde1d6,stroke:#eb6834,stroke-width:1.6px,color:#0b0b0b
    classDef plain fill:#fcfcfb,stroke:#52514e,stroke-width:1.2px,color:#0b0b0b
```

*Blue: steps · green: files you get · orange: decisions and limits.*

## What a design looks like

A scale demo: a 250-gene design on the *E. coli* K-12 genome (`design/targets/ecoli_spread_250.tsv`), followed by in-silico multiplex PCR: one simulated tube per pool, reads trimmed and aligned by `analysis/analyze.py`, shown in IGV (igv.js).

![250 targets, 250 amplicons, 4 pools](docs/igv/overview.png)

Zoom to one target: the amplicon track shows the two primers as blocks; aligned reads start inside them (primers are trimmed) and only the pool that holds this amplicon has reads.

![thrA amplicon, primers and reads](docs/igv/amplicon.png)

These reads are **simulated from the designed amplicons**, so 100% on-target and uniform coverage are true by construction. The figure shows that the design coordinates, pool assignment and analysis agree; it says nothing about real specificity or yield. The tool designs for the genes you choose; this panel is a demonstration, not a recommended panel. To rebuild the figures see `research/igv_showcase/` (unmaintained).

## Quick start

```bash
micromamba create -n amplicon_pipeline -f conda-linux-64.lock    # exact versions; or environment.yml
micromamba activate amplicon_pipeline

python design.py --project panel.example.toml                    # SARS-CoV-2, 99 overlapping tiles, ~2 min -> 2 pools
python scripts/design_viewer.py --prefix design/output/sars2_demo \
    --reference examples/sars2/NC_045512.2.fa --outdir panel_view # open panel_view/specificity.html

# after sequencing one pool:
python analysis/analyze.py --reads "raw_data/S_R{1,2}.fastq.gz" --genome examples/sars2/NC_045512.2.fa \
    --bed design/output/sars2_demo_pool_1.bed \
    --primers design/output/sars2_demo_pool_1_trimming.fasta --outdir results
python design.py --project panel.example.toml --feedback results/S/S_balancing_feedback.json
```

A project file (`panel.example.toml`) holds the options; command-line flags override it. Every design run also writes `<prefix>_manifest.json`.

**Development genome:** SARS-CoV-2 (30 kb) is the development and test genome. Design, BLAST and `check_all.sh` run in seconds to minutes, and the real-read check uses the same genome. `examples/sars2/` holds the reference and the ARTIC v3 primer scheme; `design/targets/sars2_{tiles,genes}.tsv` are the target sets. E. coli (`common_refs/`) is kept for scale tests (up to 250 targets) and the plexity study. The example design is not clean: with this tool's default primer limits only 66 of the 99 tiles design, so `panel.example.toml` sets wider limits in its `[primers]` section, and it still reports 10 failed tiles and 49 tiling gaps (3,955 bp uncovered after primer trimming). `AGENTS.md` is the operating guide for AI assistants.

## Targets

One TSV row per target: either a region (`reference`, `start`, `end`, 1-based inclusive) or a raw `sequence`.

| Column | Meaning |
|---|---|
| `name` | unique target name |
| `strand` | `+` (default) or `-`: sense strand. Design runs on the sense strand; the antisense primer is the RT primer for stranded RNA |
| `junctions` | a primer must overlap one of these joins (position of the first base after the join) |
| `span` · `avoid` | the amplicon must contain a range · no primer may overlap a range |
| `amplicon_min/max` · `max_hits` | per-target size range · allowed exact hits per primer per database |

**RNA:** design on the transcript and screen against transcriptome *and* genome (`--specificity-fasta transcripts.fa genome.fa`). The report then shows genomic amplification: on a synthetic 3-exon gene, a junction-spanning primer pair has no genomic product, while the unconstrained pair predicts a same-size one (`python scripts/make_rna_example.py`).

Ready-made sets live in `design/targets/`: `sars2_tiles.tsv` (99 overlapping tiles) and `sars2_genes.tsv` (10 ORFs, checked for start and stop codons), plus E. coli `ecoli_spread_{50,100,250}.tsv` (genome-wide, evenly spaced) and `ecoli_tiled_100.tsv`. `scripts/make_ecoli_targets.py` regenerates them.

## Design

```mermaid
%%{init: {"theme":"base","themeVariables":{"fontFamily":"Inter, system-ui, -apple-system, Segoe UI, sans-serif","fontSize":"15px","background":"#fcfcfb","primaryColor":"#cde2fb","primaryBorderColor":"#2a78d6","primaryTextColor":"#0b0b0b","secondaryColor":"#d3f1e5","tertiaryColor":"#fde1d6","lineColor":"#52514e","textColor":"#0b0b0b","clusterBkg":"#fcfcfb","clusterBorder":"#e4e3df","edgeLabelBackground":"#fcfcfb"},"flowchart":{"curve":"basis","nodeSpacing":30,"rankSpacing":42,"padding":12}}}%%
flowchart LR
    S["target +<br/>constraints"]:::step --> P["primer3<br/>candidates"]:::step
    P --> B["BLAST per<br/>database"]:::step
    B --> X["dimer and<br/>product filters"]:::step
    X --> L(["layout"]):::warn
    L -- "sparse" --> SP["fewest pools<br/>first-fit search"]:::step
    L -- "tiled" --> TP["≤ 2 pools<br/>A,C | B,D"]:::step
    SP --> O["pools · order.csv"]:::out
    TP --> O
    classDef step fill:#cde2fb,stroke:#2a78d6,stroke-width:1.6px,color:#0b0b0b
    classDef out fill:#d3f1e5,stroke:#1baf7a,stroke-width:1.6px,color:#0b0b0b
    classDef warn fill:#fde1d6,stroke:#eb6834,stroke-width:1.6px,color:#0b0b0b
    classDef plain fill:#fcfcfb,stroke:#52514e,stroke-width:1.2px,color:#0b0b0b
```

| | sparse (interspersed) | tiled (tandem) |
|---|---|---|
| Pools | fewest compatible (`--max-pools N` caps it) | **at most 2**, alternating by position |
| Cross-target product block | facing primers closer than 1000 bp | closer than the max amplicon size |
| Fails when | targets are left over under `--max-pools` | a target has no compatible pair in its pool |

`--layout auto` (default) picks tiled when target windows overlap. Tiled runs also report bases left uncovered after primer trimming (neighbours overlapping by less than two primer lengths): `<prefix>_tiling_gaps.tsv`.

**Thermodynamics:** KAPA HiFi conditions (2.5 mM Mg²⁺, 1.2 mM dNTP) on every primer3 call; primer 19-21 nt, 40-60% GC, Tm 59-61 °C by default (`--primer-min-size/-max-size`, `--primer-min-gc/-max-gc`, `--primer-min-tm/-max-tm`, or a `[primers]` section; AT-rich templates need wider limits, see Verified); dimer limits self −9, heterodimer −8, 3′ end −6 kcal/mol (`--max-*-dg`).

**Outputs** (`--output-prefix`): `_order.csv` (pool, name, sequence, role), `_pool_N.csv` (per pair, with hit counts and `predicted_products`), `_pool_N.bed`, `_pool_N_trimming.fasta`, `_manifest.json`, and `_failed_targets.txt` when relevant.

**Key flags:** `--layout`, `--max-pools`, `--max-hits`, `--amplicon-min/max`, `--primer-*` limits, `--max-cross-product`, `--fwd-tail/--rev-tail/--no-tails`, `--seed` (same inputs + seed = identical output), `--rescue` (shift the window of unconstrained region targets).

**Repeat screen:** before primer design, `--uniqueness-screen` (default `auto`) counts every k-mer (19, or the minimum primer size if smaller) of each target across the specificity FASTA(s) and fails targets where no pair of genome-unique primer sites can fit the amplicon range, with the reason. It never rejects what the full search would accept on the human test below (0 of 79 designed regions rejected) but it is a pre-filter: it caught 12 of the 21 repeat regions that failed, the other 9 fail later in the normal search. `auto` runs it for references up to about 500 Mb (seconds); for a human-size genome it takes minutes, so use `--uniqueness-screen on`. Targets eligible for `--rescue` and runs with `--force-multiprime` skip it.

**Inspect specificity:** `scripts/design_viewer.py` writes a self-contained `specificity.html` plus IGV-loadable BED tracks (untested in IGV desktop). It lists 3′-anchored near-match priming sites (≥ 16 aligned bp, ≤ 1 mismatch) that exact-hit screening cannot see.

## Analysis

```mermaid
%%{init: {"theme":"base","themeVariables":{"fontFamily":"Inter, system-ui, -apple-system, Segoe UI, sans-serif","fontSize":"15px","background":"#fcfcfb","primaryColor":"#cde2fb","primaryBorderColor":"#2a78d6","primaryTextColor":"#0b0b0b","secondaryColor":"#d3f1e5","tertiaryColor":"#fde1d6","lineColor":"#52514e","textColor":"#0b0b0b","clusterBkg":"#fcfcfb","clusterBorder":"#e4e3df","edgeLabelBackground":"#fcfcfb"},"flowchart":{"curve":"basis","nodeSpacing":30,"rankSpacing":42,"padding":12}}}%%
flowchart LR
    R["reads"]:::plain --> TR["cutadapt<br/>anchored primers"]:::step
    TR --> AL["bwa mem<br/>samtools sort"]:::step
    GN["genome"]:::plain --> AL
    AL --> M["samtools bedcov<br/>depth · stats"]:::step
    BED["amplicon BED"]:::plain --> M
    M --> RP["balance.py"]:::step
    RP --> OUT["feedback JSON<br/>summary TSV · plot"]:::out
    classDef step fill:#cde2fb,stroke:#2a78d6,stroke-width:1.6px,color:#0b0b0b
    classDef out fill:#d3f1e5,stroke:#1baf7a,stroke-width:1.6px,color:#0b0b0b
    classDef warn fill:#fde1d6,stroke:#eb6834,stroke-width:1.6px,color:#0b0b0b
    classDef plain fill:#fcfcfb,stroke:#52514e,stroke-width:1.2px,color:#0b0b0b
```

`analysis/analyze.py` needs only `cutadapt`, `bwa`, `samtools`, `numpy` and `matplotlib`. One run per pool; `{1,2}` in `--reads` marks the read number. Reads begin at the primer, so primers are matched anchored at the 5′ end (R1 forward, R2 reverse, which is what the tailed P5/P7 primers give) and unmatched reads are discarded; if none match, the run stops with an error. If your library starts reads at either end of the amplicon (ligation or tagmentation preps), add `--either-orientation`; on the real ARTIC run below the default kept 40% of pairs and this flag 80%.

Metrics: mean target coverage, on-target base %, fold-80 (mean ÷ 20th percentile of target-base depth; n/a when that percentile is 0), and % of amplicons above 0.2× the mean. On a 21-amplicon panel these reproduce the earlier Picard-based pipeline exactly, in 6 seconds instead of 1½ minutes. Output per sample in `<outdir>/<sample>/`: sorted BAM, `target_coverage.tsv`, `<sample>_balancing_feedback.json`, `_target_summary.tsv`, `_uniformity.png`, `_summary.json` (metrics plus provenance).

There is no variant calling. Amplicon reads share start positions, so over 99% of pairs look like duplicates: do not run duplicate-aware tools or filters on the BAM.

## Feedback loop

```mermaid
%%{init: {"theme":"base","themeVariables":{"fontFamily":"Inter, system-ui, -apple-system, Segoe UI, sans-serif","fontSize":"15px","background":"#fcfcfb","primaryColor":"#cde2fb","primaryBorderColor":"#2a78d6","primaryTextColor":"#0b0b0b","secondaryColor":"#d3f1e5","tertiaryColor":"#fde1d6","lineColor":"#52514e","textColor":"#0b0b0b","clusterBkg":"#fcfcfb","clusterBorder":"#e4e3df","edgeLabelBackground":"#fcfcfb"},"flowchart":{"curve":"basis","nodeSpacing":30,"rankSpacing":42,"padding":12}}}%%
flowchart LR
    C["per-target coverage"]:::plain --> R["R = target / panel mean"]:::step
    R --> Z(["R"]):::warn
    Z -- "0" --> RD["redesign<br/>(--rescue)"]:::warn
    Z -- "< 0.5" --> UP["raise ×1/R<br/>(max 10)"]:::step
    Z -- "0.5 – 2" --> KP["keep"]:::out
    Z -- "> 2" --> DN["lower ×1/R<br/>(min 0.1)"]:::step
    UP --> REC["rebalancing_recipe.csv<br/>10 µL × factor"]:::out
    DN --> REC
    classDef step fill:#cde2fb,stroke:#2a78d6,stroke-width:1.6px,color:#0b0b0b
    classDef out fill:#d3f1e5,stroke:#1baf7a,stroke-width:1.6px,color:#0b0b0b
    classDef warn fill:#fde1d6,stroke:#eb6834,stroke-width:1.6px,color:#0b0b0b
    classDef plain fill:#fcfcfb,stroke:#52514e,stroke-width:1.2px,color:#0b0b0b
```

## Reproducibility

* **`conda-linux-64.lock`** pins 200 packages to exact builds (primer3-py 2.3.1, BLAST 2.17.0, bwa 0.7.19, samtools 1.25, cutadapt 5.2). `environment.yml` is the loose spec for other platforms.
* **Design manifest:** input hashes, every option, the seed, tool versions, git commit and a dirty-tree flag. **Analysis `_summary.json`:** the command, tool versions and input hashes.
* The same inputs and seed give identical output (three real runs of a 92-target panel were byte-identical).
* **Honest gaps:** the lock file is Linux x86-64 only and there is no container. Version drift rarely matters but is not zero: primer3 2.3.1 and 2.4.0 differ on 3 of 1,300 ΔG values (up to 0.9 kcal/mol, enough to flip a pair near the −6 limit), though they designed identical primers. Real reads are not in the repository; simulated ones are reproducible from `--seed`.
* Workflow engines are not needed for this: the previous Nextflow pipeline pinned no tools either.

## Verified

Real BLAST+, bwa, samtools and cutadapt, simulated reads. Run everything with `bash scripts/check_all.sh` (about 2 minutes); CI runs lint and unit tests on every push, and `real-tools.yml` runs the full check on demand.

| Test | Result |
|---|---|
| SARS-CoV-2 ORFs, sparse (`sars2_genes`), wider primer limits | 10/10 designed, 1 pool |
| **SARS-CoV-2 tiled, 99 tiles, tool-default primer limits** | **66/99 designed**, 33 failed, **11,680 bp (39%) uncovered after primer trimming**. The defaults (19-21 nt, GC 40-60%) cannot place primers on a 38% GC genome; the 196 real ARTIC primers are 22-30 nt |
| SARS-CoV-2 tiled, same tiles, wider limits (`panel.example.toml`: 20-30 nt, GC >= 30%, Tm 57-64) | 89/99 designed in 2 pools [45, 44], 10 failed, 3,955 bp uncovered (49 stretches), ~2.4 min. ARTIC covers the genome with 98 amplicons. Neighbour overlap is not optimised |
| Do our metrics predict real ARTIC amplicon coverage? (`research/benchmark_sars2/`) | No. Real coverage spans 44x-1522x; Spearman correlation with 12 metrics (Tm, GC, dimers, 3' stability, amplicon length and GC, ...) is |rho| <= 0.19, none significant at p < 0.05 (one dataset, one chemistry) |
| 100 genes spread over the genome (`ecoli_spread_100`) | 100/100 designed, 3 pools [58, 37, 5], 128 s |
| 250 genes spread over the genome (`ecoli_spread_250`) | 250/250 designed, 4 pools [89, 83, 53, 25], no failures, no unintended predicted products; about 11 min with 5 packing iterations. Simulated reads: every pool mean about 380×, fold-80 1.26–1.29 (see IGV above) |
| Tiled design, 100 contiguous tiles (`ecoli_tiled_100`) | **Only 75/100 designed** (25 failed: no compatible pair in the tiled pool), 41 tiling gaps. A dense tiled panel is hard for this packer; read the failure list before ordering |
| RNA junction example | junction-spanning pair: no genomic product; unconstrained pair: same-size genomic product |
| Analysis, 21-gene panel | 6 s; 6,299 of 6,300 pairs kept; mean 397×, fold-80 1.33; identical to the Picard run |
| **Real reads:** public MiSeq ARTIC v3 SARS-CoV-2 run (98 amplicons, 2 tiled pools, 27.7k pairs); `bash scripts/real_data_check.sh`, about 1 min, needs internet | 80% of pairs kept with `--either-orientation` (40% without), 99.96% mapped, mean 409×, **fold-80 5.6** (simulated: 1.3); all 98 amplicons covered, lowest 44× (amp98), with the known weak ARTIC v3 amplicons 64 and 70 among the lowest five; the recipe asks for 0.27× to 9× concentration changes |
| Human GRCh38 (primary assembly), 100 random 400 bp regions with simple-repeat and GC filtering; **not a panel, arbitrary regions** | 79 designed in 3 pools [45, 26, 8], 18 min; the 21 failed regions are high-copy repeats (median about 26,000 genomic matches vs 1). The repeat screen alone took 6.4 min on this genome and caught 12 of the 21 |
| Unit tests | 84 (82 plus two end-to-end tests that skip without bwa/samtools) |

**Limits:** the simulated-read results are uniform by construction. The analysis has been run on one real public dataset, which is viral, single-sample and not designed by this tool, so the rebalancing recipe's effect is untested, and no human amplicon reads have been run (the hosts that hold them are unreachable from the dev sandbox); targets shorter than the minimum amplicon cannot be designed; BLAST screening counts exact hits (near-matches are informational); packing is a heuristic; the 1000 bp cross-product distance is a chosen default, not a validated one.

## In-silico plexity study

How far does a tiled design go as plexity grows, and what does chemistry buy back? **935 overlapping E. coli tiles** (of 1,000; 65 had no specific primer pair) in two alternating pools, real primer3 and BLAST, 10 to 935 targets. Figures 1, 2 and 4 carry the findings; 3 and 5 are illustrative (model-based and E. coli-only). Code and details: [`research/plexity_study/`](research/plexity_study/README.md) (unmaintained).

![Primer-dimer potential vs plexity](research/plexity_study/figures/fig1_dimer_potential_vs_plexity.png)

![Capacity vs plexity](research/plexity_study/figures/fig2_capacity_vs_plexity.png)

![Chemistry ladder at 935-plex](research/plexity_study/figures/fig4_rescue_ladder.png)

![Uniformity vs plexity](research/plexity_study/figures/fig3_uniformity_vs_plexity.png)

![Off-target vs plexity](research/plexity_study/figures/fig5_offtarget_vs_plexity.png)

| | 100-plex | 935-plex |
|---|---|---|
| Targets placeable under the plain dimer limit (first-fit / + repair) | 75% / 77% | 26% / 31% |
| Risky dimer pairs: plain / + hot-start / full ladder | 155 / 2 / 0 | 15,750 / 205 / 3 |
| Mean primer sequestered by dimers: plain / + hot-start / full ladder | 10.1% / 0.07% / 0.004% | 62% / 0.85% / 0.04% |
| Competitive off-target sites: plain / + hot-start / full ladder | 367 / 100 / 9 | 3,495 / 1,096 / 126 |

Brute-force screening alone places roughly 145 targets per pool under the plain limit; the levers that rescue it are commodity ones (hot-start polymerase, higher annealing temperature, lower primer concentration, rebalanced concentrations), though I cannot vouch for their IP status. **Read it as a model:**

1. The plain baseline holds 250 nM per primer at every plexity (0.23 mM total at 935-plex), a worst case; most commercial mixes are already hot-start.
2. Hot-start is modelled as dimers extending at 55 °C instead of 37 °C, the single assumption behind the biggest effect; real high-plex PCR still makes dimers with hot-start.
3. "Risky" means pair occupancy ≥ 0.42%, calibrated to the tool's own −6 kcal/mol limit, not to data; no sensitivity analysis has been run.
4. The uniformity model is uncalibrated; off-target covers a 4.6 Mb genome (a human-scale genome has ~1,000× more candidate sites).
5. Greedy selection plus one-swap repair is a pessimistic bound on what brute force can place.

## Repository

```
design.py                    design, packing, feedback, project file, manifest
core/biophysics.py           primer3 thermodynamics (KAPA conditions)
analysis/                    analyze.py, metrics.py, balance.py, plot_uniformity.py, simulator
core/uniqueness.py          genome k-mer copy counts for the repeat screen
scripts/                     design_viewer.py, check_all.sh, real_data_check.sh, make_ecoli_targets.py, ...
examples/sars2/              SARS-CoV-2 reference + ARTIC v3 scheme (the development genome)
design/targets/              sars2_* and ecoli_* target sets      common_refs/   E. coli K-12 genome + GFF (scale tests)
research/                    plexity_study/, igv_showcase/, benchmark_sars2/ (unmaintained)
docs/igv/                    IGV figures
test_design.py               unit tests               conda-linux-64.lock · environment.yml
```

Earlier iterations removed (recoverable from git history): the tail-only mode, Docker, Nextflow, GATK variant calling, Picard, FastQC, MultiQC, the hairpin clamp, and the HTML dashboard.
