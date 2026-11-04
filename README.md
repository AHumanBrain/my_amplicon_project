# 🧬 Multiplex Amplicon Toolkit

Design multiplex PCR primer panels, simulate reads, analyse sequencing data with Nextflow, and feed coverage back into pool rebalancing.

| | |
|---|---|
| **Design** | `design.py` → primer3 + BLAST+ + dimer screening + pool packing |
| **Analyse** | `analysis/main.nf` → Cutadapt, BWA, Picard, GATK, deepTools, MultiQC |
| **Close the loop** | coverage JSON → pipetting recipe + redesign of failed targets |
| **Stack** | conda env (`environment.yml`) or Docker (`Dockerfile`); Python 3.11, Nextflow ≥ 25.04 |

---

## 1. Big picture

```mermaid
flowchart LR
    subgraph IN["Inputs"]
        G["Genome FASTA"]
        A["GFF annotation"]
        T["Target gene list"]
    end
    subgraph DES["design.py"]
        D1["primer3 candidates"] --> D2["BLAST specificity"] --> D3["Dimer filters"] --> D4["Pool packing"]
    end
    subgraph OUT["Design outputs"]
        O1["*.csv order sheet"]
        O2["*.bed amplicons"]
        O3["*_trimming.fasta"]
    end
    subgraph SEQ["Data"]
        W["Wet-lab FASTQ"]
        S["Simulated FASTQ<br/>generate_simulated_fastq.py"]
    end
    subgraph NF["analysis/main.nf"]
        N1["Trim → Align → Metrics → Variants"]
    end
    subgraph FB["Feedback"]
        F1["balancing_feedback.json"]
        F2["rebalancing_recipe.csv"]
    end

    G & A & T --> D1
    D4 --> O1 & O2 & O3
    O2 --> S
    S --> N1
    W --> N1
    O2 & O3 --> N1
    N1 --> F1 --> F2
    F1 -- "REDESIGN flags" --> D1
```

---

## 2. Quick start

```bash
# environment (conda/micromamba) – or: docker build -t amplicon-project .
micromamba create -y -n amplicon_pipeline -f environment.yml
micromamba activate amplicon_pipeline

# 1. design a panel (21 housekeeping genes ≈ 7 s)
python design.py --genome common_refs/ecoli_genome.cleaned.fna \
    --gff common_refs/genomic.gff \
    --target-file design/targets/housekeeping_genes.txt \
    --blast-db ncbi_data/ecoli_db --output-prefix design/output/hk

# 2. simulate reads (also writes raw_data/HK.truth.tsv, the injected SNVs)
python analysis/generate_simulated_fastq.py --genome common_refs/ecoli_genome.cleaned.fna \
    --bed design/output/hk_pool_1.bed --output-dir raw_data --sample-name HK --coverage 300 --seed 1

# 3. analyse
nextflow run analysis/main.nf --genome common_refs/ecoli_genome.cleaned.fna \
    --reads "raw_data/HK_R{1,2}.fastq.gz" --bed design/output/hk_pool_1.bed \
    --primers design/output/hk_pool_1_trimming.fasta --outdir results

# 4. rebalance / redesign from coverage
python design.py --feedback results/coverage/HK_balancing_feedback.json --rescue \
    --output-prefix design/output/hk_v2
```

---

## 3. Module 1 – Primer design (`design.py`)

### 3.1 Per-target logic

```mermaid
flowchart TD
    S["Extract gene sequence from GFF + FASTA"] --> P["primer3: N candidate pairs<br/>Strategy 0 ideal"]
    P --> B["ONE blastn call per candidate batch<br/>exact full-length hits only"]
    B --> C1["BLAST: at most 1 genomic hit<br/>fails closed on BLAST error"]
    C1 --> C2["Self-dimer ≥ -9.0 kcal/mol"]
    C2 --> C3["F/R heterodimer ≥ -8.0"]
    C3 --> C4["3' end stability ≥ -6.0 both directions"]
    C4 --> OK["Keep pair"]
    P -. "no valid pair" .-> FB1["Fallback strategies 1 to 3<br/>longer / shorter / relaxed Tm"]
    FB1 -. "still none" .-> R["--rescue: shift window ±25 bp, ±45 bp"]
    R -. "still none" .-> X["Target reported as failed"]
```

### 3.2 Packing candidates into pools

```mermaid
flowchart LR
    M["Candidates per target<br/>sorted by strategy, then numeric rank"] --> I["Monte Carlo first-fit<br/>iteration 0 = genomic order,<br/>others shuffled with --seed"]
    I --> Q{"Fits an existing pool?<br/>cross-dimer + overlap checks"}
    Q -- yes --> A["Add best-ranked compatible candidate"]
    Q -- no --> N["Open new pool"]
    A & N --> K["Keep run with fewest pools"]
    K --> W["Write *_pool_N.csv / .bed / _trimming.fasta / .dimers.txt"]
```

### 3.3 Thermodynamic model

| Setting | Value | Note |
|---|---|---|
| Mg²⁺ / dNTP / Na⁺ / DNA | 2.5 mM / 1.2 mM / 50 mM / 50 nM | KAPA HiFi; **now applied to every primer3 call** |
| Tm window | 59–61 °C (fallback 55–65) | GC 40–60 %, 19–21 nt |
| Amplicon | 150–250 bp (fallback 100–350) | must fit inside the gene |
| Self-dimer / heterodimer / 3′ end | −9.0 / −8.0 / −6.0 kcal/mol | override with `--max-*-dg` |

### 3.4 Modes and key flags

```mermaid
flowchart TD
    A{"Arguments"} -- "--feedback" --> M3["Mode 3: rebalance + redesign<br/>(takes priority)"]
    A -- "--target-file" --> M1["Mode 1: full design"]
    A -- "--tail-fwd-file + --tail-rev-file" --> M2["Mode 2: add adapter tails only"]
```

| Flag | Default | Purpose |
|---|---|---|
| `--num-candidates` | 25 | pairs designed per target |
| `--monte-carlo-iters` | 20 | packing restarts |
| `--seed` | 0 | reproducible packing |
| `--force-single-pool` | off | one tube; targets that don't fit are **listed in the log**, not silently dropped |
| `--force-sparse` | off | enforce zero amplicon overlap |
| `--force-multiprime` | off | allow multi-copy loci (e.g. 16S) |
| `--oligo-format` | `rc_tailed` | or `fwd_tailed` (+ `--add-hairpin-clamp`) |
| `--output-prefix` | `design/output/final_primers` | |
| `--log-file` | none | optional log file |

---

## 4. Module 2 – Nextflow pipeline (`analysis/main.nf`)

```mermaid
flowchart TD
    R["Reads R1/R2"] --> FQ["FASTQC"]
    R --> TR["TRIM_PRIMERS<br/>cutadapt: R1 fwd primers, R2 rev primers"]
    GN["Genome"] --> BI["BUILD_INDICES<br/>bwa, faidx, dict"]
    BED["Amplicon BED"] --> PI["PREP_INTERVALS"]
    BI --> PI
    TR --> AL["ALIGN<br/>bwa mem piped to samtools sort"]
    BI --> AL
    AL --> SI["SORT_INDEX_BAM<br/>index only"]
    SI --> MD["MARK_DUPLICATES<br/>informational"]
    SI --> PM["PICARD_METRICS"]
    PI --> PM
    SI --> VC["GERMLINE_CALLING or MUTECT2<br/>--mode germline / somatic"]
    SI --> BW["COVERAGE_TRACKS .bw"]
    PM --> PC["PLOT_COVERAGE<br/>HTML + feedback JSON + uniformity PNG"]
    FQ & TR & PM & PC & MD --> MQ["MULTIQC"]
```

### Primer trimming

```mermaid
flowchart LR
    PF["*_trimming.fasta<br/>tail+primer, names end _F / _R"] --> SPLIT["awk split"]
    SPLIT --> FW["fwd_primers.fa"]
    SPLIT --> RV["rev_primers.fa"]
    FW -- "-g" --> CA["cutadapt -O 15<br/>--discard-untrimmed"]
    RV -- "-G" --> CA
    R1["R1"] --> CA
    R2["R2"] --> CA
    CA --> TRIMMED["trimmed pairs"]
```

Reads start inside tail+primer, so partial 5′ matches are expected; `--min_primer_overlap` (default 15) prevents spurious short hits.

### Why variants are called on the *non*-duplicate-marked BAM

```mermaid
flowchart LR
    A["Amplicon reads share start positions"] --> B["Picard: 99.7 % duplicates"]
    B --> C{"BAM given to GATK"}
    C -- "duplicate-marked" --> D["GATK drops dups: mean DP 2,<br/>17 variants (wrong)"]
    C -- "sorted, unmarked" --> E["mean DP ~83 at calling sites:<br/>27 of 42 SNVs, 0 false positives"]
```

### Outputs (`--outdir`)

| Path | Content |
|---|---|
| `trimmed_reads/` | trimmed FASTQ + cutadapt log |
| `bams/` | sorted + indexed BAM, marked BAM, dup metrics |
| `qc/picard/`, `qc/fastqc/` | Picard PCR/alignment/insert metrics, FastQC |
| `variants/germline` or `somatic/` | GVCF / Mutect2 VCF |
| `coverage/` | bigWig, `*_coverage_report.html`, `*_uniformity.png`, `*_balancing_feedback.json` |
| `multiqc/` | aggregated report |

---

## 5. Module 3 – Feedback loop

```mermaid
flowchart LR
    COV["Per-target mean coverage"] --> REL["R = target / panel mean"]
    REL --> ACT{"Action"}
    ACT -- "R = 0" --> RD["REDESIGN<br/>gene ID recovered from amplicon name,<br/>re-run design with --rescue"]
    ACT -- "R below 0.5" --> UP["INCREASE ×1/R (max 10)"]
    ACT -- "0.5 to 2" --> KEEP["Maintain"]
    ACT -- "R above 2" --> DOWN["DECREASE ×1/R (min 0.1)"]
    UP & DOWN & KEEP --> REC["rebalancing_recipe.csv<br/>10 µL × factor"]
```

---

## 6. Simulator (`analysis/generate_simulated_fastq.py`)

Deterministic for a given `--seed`; streams to gzip. Options: `--coverage`, `--non-uniform` (10–200 % per amplicon), `--thermal-bias` (coverage falls 20 %/°C from a per-target optimum near 62 °C, stable via CRC32), `--mutation-rate` (default 1 %), `--error-rate` (0.1 %). Writes `<sample>.truth.tsv` listing every injected SNV.

---

## 7. What changed (this revision)

```mermaid
flowchart TD
    subgraph BROKEN["Found broken"]
        B1["design.py never imported primer3 / csv;<br/>adapter + hairpin constants undefined"]
        B2["Bare excepts hid errors as 'no dimer'"]
        B3["Relative paths as Nextflow path inputs:<br/>CI stub run failed"]
        B4["GATK ran on duplicate-marked amplicon BAM"]
        B5["Picard '?' metric crashed the report"]
        B6["Salt calibration documented but never applied"]
    end
    subgraph FIXED["Fixed"]
        F1["Imports + constants restored, pyflakes in CI"]
        F2["Specific exceptions, invalid seqs rejected, BLAST fails closed"]
        F3["Reference inputs via Channel.fromPath"]
        F4["Variant calling on sorted BAM"]
        F5["to_float() helper"]
        F6["KAPA conditions passed to primer3"]
    end
    B1 --> F1
    B2 --> F2
    B3 --> F3
    B4 --> F4
    B5 --> F5
    B6 --> F6
```

**Also improved:** BLAST batched (1 call per batch instead of 2 per pair); numeric candidate ranking (was a string sort); Picard-safe BED names; thresholds propagated to worker processes; `--seed`; single-pool no longer drops targets silently; `bwa | samtools sort` with `task.cpus`; correct R1/R2 primer trimming; no `errorStrategy 'ignore'`; seeded simulator with truth file; numpy k-mer script (verified against brute force; 4.6 Mbp, k=12–16 in ~7 s); HTML-escaped reports; `.gitignore` no longer contradicts tracked files (example panels in `design/examples/`).

---

## 8. Verified results

Run with real BLAST+, BWA, Picard, GATK and Cutadapt from `environment.yml`.

| Test | Result |
|---|---|
| 21-gene panel | designed in ~7 s → **1 pool**; all 42 primers have exactly 1 perfect genomic hit |
| 100-gene stress list | 92 designed in **3 pools (62 / 27 / 3)**; 8 failed (e.g. `thrL`, `sokC` are shorter than the minimum amplicon) |
| Primer trimming | 6,289 / 6,300 pairs retained (99.8 %) |
| Pipeline | all 12 processes succeed on panel reads and on the CI fixture |
| Panel coverage | mean 390×, fold-80 1.32, 100 % of targets above 0.2× |
| Unit tests | 27 passing |

```mermaid
pie showData
    title "42 injected SNVs: outcome"
    "Called by GATK" : 27
    "Missed: under a primer (trimmed by design)" : 10
    "Missed: amplicon interior (unexplained)" : 5
```

---

## 9. Known limitations / open items

* **5 interior SNVs are not called** although the BAM pileup shows them at 100 % allele fraction and 300×; under investigation.
* **Short genes** (below ~100 bp) cannot fit an amplicon; designing on gene ± flank would fix this.
* **Specificity** counts only exact full-length matches; 3′-mismatch near-hits are not scored.
* **Packing** is randomised first-fit; graph colouring or ILP might beat 3 pools on the 100-gene list.
* **Assumption:** reads begin at the primer (adapter tail not sequenced). Revisit `--min_primer_overlap` if your library prep differs.

---

## 10. Repository layout and CI

```mermaid
flowchart TD
    ROOT["repo"] --> DP["design.py"]
    ROOT --> CORE["core/<br/>biophysics · adapters · io_utils"]
    ROOT --> AN["analysis/<br/>main.nf · nextflow.config · plot_*.py · generate_simulated_fastq.py · analysis_utils.py"]
    ROOT --> SC["scripts/<br/>curate_sparse_targets · measure_kmer_uniqueness · make_test_fixture"]
    ROOT --> DS["design/<br/>targets/ · examples/ · output/ (ignored)"]
    ROOT --> CR["common_refs/<br/>E. coli genome + GFF"]
    ROOT --> TS["test_design.py"]
    ROOT --> CI[".github/workflows/test-pipeline.yml"]
```

```mermaid
flowchart LR
    P["push / PR to main"] --> J1["Job 1: pyflakes + pytest"]
    J1 --> J2["Job 2: micromamba env"]
    J2 --> FX["make_test_fixture.py<br/>8 kb genome slice, 3 amplicons, simulated reads"]
    FX --> ST["Nextflow stub run"]
    ST --> RL["Real run: BWA / Picard / GATK"]
    RL --> CK["Assert 6 key outputs exist"]
```

---

## 11. FAQ

**Why several pool files?** Overlapping amplicons or mutually dimerising primers can't share a tube; the packer splits them (`_pool_1`, `_pool_2`, …).
**Fold-80 penalty?** Extra sequencing needed to bring 80 % of targets to the mean; lower is more uniform.
**Single-end data?** Detected automatically; Cutadapt/BWA adapt, insert-size metrics are skipped gracefully.
**`--rescue`?** Shifts the extraction window (±25 bp, ±45 bp) when no pair passes at the default position.
