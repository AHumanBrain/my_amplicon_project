#!/bin/bash

# run_analysis.sh: Wrapper to run the Amplicon Analyzer pipeline
# Points to shared references and design outputs in the unified structure.

# Find the genome file dynamically
GENOME_FILES=(../common_refs/*.cleaned.fna)
if [[ ${#GENOME_FILES[@]} -eq 0 ]] || [[ ! -e "${GENOME_FILES[0]}" ]]; then
    echo "Error: No genome file found matching ../common_refs/*.cleaned.fna"
    exit 1
elif [[ ${#GENOME_FILES[@]} -gt 1 ]]; then
    echo "Error: Multiple genome files found in ../common_refs/: ${GENOME_FILES[*]}. Please ensure only one cleaned genome exists."
    exit 1
fi
GENOME=$(realpath "${GENOME_FILES[0]}")

# Set other variables relative to the analysis directory
# Note: Pointing to pool_1 by default.
PRIMERS=$(realpath "../design/output/final_primers_pool_1_trimming.fasta")
BED=$(realpath "../design/output/final_primers_pool_1.bed")
READS=$(realpath "../raw_data")/*_R{1,2}.fastq.gz
OUTDIR=$(realpath -m "../results") # -m because it might not exist yet

if [[ ! -f "$PRIMERS" ]]; then 
    echo "Warning: Primers file not found at $PRIMERS. Ensure you have run the design script."
fi

if [[ ! -f "$BED" ]]; then 
    echo "Warning: BED file not found at $BED. Ensure you have run the design script."
fi

# Run Nextflow
# Note: --reads uses a glob pattern that bash might expand, so we pass it in quotes.
nextflow run main.nf \
    -profile conda \
    --genome "$GENOME" \
    --primers "$PRIMERS" \
    --bed "$BED" \
    --reads "$READS" \
    --outdir "$OUTDIR" \
    -resume
