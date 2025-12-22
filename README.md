# My Amplicon Project (Unified Design & Analysis)

This project integrates the **Multiplex PCR Primer Designer** and the **Amplicon Analyzer** pipeline into a single, cohesive workflow. It automates the entire process from selecting gene targets to calling variants (VCF) from sequencing data.

---

## 📂 Project Structure

```plaintext
my_amplicon_project/
├── common_refs/               <-- SHARED TRUTH (Reference Genomes)
│   ├── ecoli_genome.cleaned.fna (Unix line endings required)
│   └── genomic.gff
├── design/                    <-- PRIMER DESIGN TOOL
│   ├── design_v7.9.py         (Python script)
│   ├── requirements.txt
│   ├── targets/               (Input gene lists)
│   └── output/                (Generated CSVs, BEDs, Fastas)
├── analysis/                  <-- NEXTFLOW PIPELINE
│   ├── setup_wsl.sh           (One-click Enviroment Setup)
│   ├── main.nf                (Nextflow Pipeline)
│   ├── nextflow.config
│   ├── generate_simulated_fastq.py
│   └── run_analysis.sh        (Execution Wrapper)
├── raw_data/                  <-- INPUT FASTQ FILES
├── results/                   <-- FINAL OUTPUTS (VCF, HTML Reports)
└── ncbi_data/                 <-- BLAST DATABASE
```

---

## 🚀 Getting Started (Windows/WSL)

Since this project uses high-performance bioinformatics tools (BWA, Samtools, GATK), it runs inside **WSL (Windows Subsystem for Linux)**.

### 1. Environment Setup (One-Time)
We have provided a script to automatically install **Miniconda** and all required dependencies (Nextflow, GATK, etc.).

1.  Open your WSL terminal (e.g., Ubuntu).
2.  Navigate to the project directory:
    ```bash
    cd /mnt/c/Users/YOUR_USER/Documents/my_amplicon_project
    ```
3.  Run the setup script:
    ```bash
    wsl bash analysis/setup_wsl.sh
    ```
    *This will install Miniconda and create the `amplicon_pipeline` environment.*

---

## 🧬 workflow

### Step 1: Design Primers
Run the Python design script to generate primers for your targets. This runs on Windows (or WSL) and uses the shared reference genome.

```bash
cd design
python design_v7.9.py --target-file targets/housekeeping_genes.txt --oligo-format fwd_tailed --add-hairpin-clamp
```
*   **Outputs** (in `design/output/`):
    *   `final_primers.csv`: Order sheet.
    *   `final_primers.bed`: Targets for analysis.
    *   `final_primers_trimming.fasta`: Primer sequences for trimming.

### Step 2: Acquire Data (Real or Simulated)
Place your paired-end FASTQ files in `raw_data/` (e.g., `Sample_R1.fastq.gz`, `Sample_R2.fastq.gz`).

**Option B: Generate Simulated Test Data**
If you don't have real data yet, generate synthetic reads with mutations:
```bash
wsl python analysis/generate_simulated_fastq.py
```
*   This creates `Simulated_Sample_R{1,2}.fastq.gz` in `raw_data/`.

### Step 3: Run Analysis
Execute the Nextflow pipeline inside WSL. The wrapper script handles path configuration.

```bash
# 1. Activate Environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate amplicon_pipeline

# 2. Go to Analysis folder
cd analysis

# 3. Run Pipeline
bash run_analysis.sh
```

---

## � Outputs (`results/`)

*   **`variants/*.vcf.gz`**: Final called variants (VCF format).
*   **`multiqc/multiqc_report.html`**: Interactive quality control report.
*   **`bams/*.sorted.bam`**: Aligned reads (viewable in IGV).
*   **`qc/`**: FastQC and Picard metrics.

---

## 🛠️ Troubleshooting

*   **"Bad input / Non-standard base" (GATK Error)**:
    *   Cause: Reference genome has Windows line endings (`CRLF`).
    *   Fix: Run `dos2unix common_refs/ecoli_genome.cleaned.fna` or `tr -d '\r' < original > fixed`.
*   **"Process requirement exceeds available memory"**:
    *   Cause: WSL has limited RAM by default.
    *   Fix: Edit `analysis/nextflow.config` and reduce memory to `'4 GB'`.
