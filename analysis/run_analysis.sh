#!/bin/bash

# run_analysis.sh: Wrapper to run the Amplicon Analyzer pipeline
# Supports config-based execution via params.json

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Check for params file argument or use default
PARAMS_FILE="${1:-params.json}"

if [[ -f "$PARAMS_FILE" ]]; then
    echo "Using parameters from: $PARAMS_FILE"
    nextflow run main.nf \
        -profile conda \
        -params-file "$PARAMS_FILE" \
        -resume
else
    echo "No params file found at $PARAMS_FILE. Using command-line defaults."
    
    # Find the genome file dynamically
    GENOME_FILES=(../common_refs/*.cleaned.fna)
    if [[ ${#GENOME_FILES[@]} -eq 0 ]] || [[ ! -e "${GENOME_FILES[0]}" ]]; then
        echo "Error: No genome file found matching ../common_refs/*.cleaned.fna"
        exit 1
    elif [[ ${#GENOME_FILES[@]} -gt 1 ]]; then
        echo "Error: Multiple genome files found. Please ensure only one cleaned genome exists."
        exit 1
    fi
    GENOME=$(realpath "${GENOME_FILES[0]}")

    PRIMERS=$(realpath "../design/output/final_primers_pool_1_trimming.fasta")
    BED=$(realpath "../design/output/final_primers_pool_1.bed")
    READS=$(realpath "../raw_data")/*_R{1,2}.fastq.gz
    OUTDIR=$(realpath -m "../results")

    MODE="${2:-germline}"
    echo "Running in mode: $MODE"

    nextflow run main.nf \
        -profile conda \
        --genome "$GENOME" \
        --primers "$PRIMERS" \
        --bed "$BED" \
        --reads "$READS" \
        --outdir "$OUTDIR" \
        --mode "$MODE" \
        -resume
fi
