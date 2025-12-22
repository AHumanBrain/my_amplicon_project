#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// --- Input Parameters ---
params.reads = "$baseDir/../raw_data/*_R{1,2}.fastq.gz"
params.genome = "$baseDir/../common_refs/*.cleaned.fna"
params.primers = "$baseDir/../design/output/final_primers_pool_1_trimming.fasta" 
params.bed = "$baseDir/../design/output/final_primers_pool_1.bed"
params.outdir = "$baseDir/../results"

log.info """
A M P L I C O N   S E Q   P I P E L I N E
=========================================
genome       : ${params.genome}
reads        : ${params.reads}
primers      : ${params.primers}
targets (BED): ${params.bed}
outdir       : ${params.outdir}
"""

// --- Processes ---

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.html", emit: html
    path "*.zip", emit: zip

    script:
    """
    fastqc -q $reads
    """
}

process TRIM_PRIMERS {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path primers

    output:
    tuple val(sample_id), path("*.trimmed.fastq.gz"), emit: trimmed_reads
    path "*.log", emit: log

    // Using -g (5' adapter) for R1 and -G (5' adapter) for R2
    // Assuming primers are anchored at the 5' end.
    script:
    """
    cutadapt \
    -g file:${primers} \
    -G file:${primers} \
    -o ${sample_id}_R1.trimmed.fastq.gz \
    -p ${sample_id}_R2.trimmed.fastq.gz \
    --discard-untrimmed \
    ${reads[0]} ${reads[1]} > ${sample_id}_cutadapt.log
    """
}

process BUILD_INDICES {
    tag "Genome Indices"
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path genome

    output:
    path "${genome}.*", emit: bwa_index
    path "${genome}.fai", emit: fai
    path "${genome.baseName}.dict", emit: dict
    // We emit the genome file again to pass it along with indices
    path genome, emit: fasta 

    script:
    """
    bwa index $genome
    samtools faidx $genome
    picard CreateSequenceDictionary R=$genome O=${genome.baseName}.dict
    """
}

process ALIGN {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)
    path index_files // dummy input to ensure index is ready
    path genome

    output:
    tuple val(sample_id), path("*.bam"), emit: bam

    script:
    """
    bwa mem -t 4 -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA' $genome ${reads[0]} ${reads[1]} | \
    samtools view -b - > ${sample_id}.bam
    """
}

process SORT_INDEX_BAM {
    tag "$sample_id"
    publishDir "${params.outdir}/bams", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: sorted_bam

    script:
    """
    samtools sort -o ${sample_id}.sorted.bam $bam
    samtools index ${sample_id}.sorted.bam
    """
}

process PREP_INTERVALS {
    // Converts simple BED file to Picard IntervalList
    input:
    path bed
    path dict

    output:
    path "targets.interval_list", emit: interval_list

    script:
    """
    picard BedToIntervalList \
    I=$bed \
    O=targets.interval_list \
    SD=$dict
    """
}

process PICARD_METRICS {
    tag "$sample_id"
    publishDir "${params.outdir}/qc/picard", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path interval_list

    output:
    path "*_metrics.txt"

    script:
    """
    # 1. Alignment Metrics
    picard CollectAlignmentSummaryMetrics \
    R=$fasta I=$bam O=${sample_id}_alignment_metrics.txt

    # 2. Insert Size Metrics
    picard CollectInsertSizeMetrics \
    I=$bam O=${sample_id}_insert_size_metrics.txt H=${sample_id}_insert_size_histogram.pdf

    # 3. Targeted PCR Metrics
    picard CollectTargetedPcrMetrics \
    I=$bam O=${sample_id}_pcr_metrics.txt \
    AMPLICON_INTERVALS=$interval_list \
    TARGET_INTERVALS=$interval_list \
    R=$fasta
    """
}

process VARIANT_CALLING {
    tag "$sample_id"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path fai
    path dict

    output:
    path "${sample_id}.vcf.gz", emit: vcf

    // Using PCR indel model conservatively
    script:
    """
    gatk HaplotypeCaller \
    -R $fasta \
    -I $bam \
    -O ${sample_id}.vcf.gz \
    --pcr-indel-model CONSERVATIVE \
    -ERC GVCF
    """
}

process COVERAGE_TRACKS {
    tag "$sample_id"
    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "*.bw"

    script:
    """
    bamCoverage -b $bam -o ${sample_id}.bw
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path '*' // Collects all files from previous processes

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

// --- Workflow ---

workflow {
    // 1. Channel for reads
    Channel
        .fromFilePairs(params.reads)
        .set { read_pairs_ch }

    // 2. Build Reference Indices (runs once)
    BUILD_INDICES(params.genome)
    
    // 3. Prepare Intervals (runs once)
    PREP_INTERVALS(params.bed, BUILD_INDICES.out.dict)

    // 4. Raw QC
    FASTQC(read_pairs_ch)

    // 5. Trim Primers
    TRIM_PRIMERS(read_pairs_ch, params.primers)

    // 6. Align
    ALIGN(TRIM_PRIMERS.out.trimmed_reads, BUILD_INDICES.out.bwa_index, params.genome)

    // 7. Sort & Index
    SORT_INDEX_BAM(ALIGN.out.bam)

    // 8. Metrics & Visualization
    PICARD_METRICS(
        SORT_INDEX_BAM.out.sorted_bam, 
        BUILD_INDICES.out.fasta, 
        PREP_INTERVALS.out.interval_list
    )
    
    COVERAGE_TRACKS(SORT_INDEX_BAM.out.sorted_bam)

    // 9. Variant Calling
    VARIANT_CALLING(
        SORT_INDEX_BAM.out.sorted_bam, 
        BUILD_INDICES.out.fasta, 
        BUILD_INDICES.out.fai, 
        BUILD_INDICES.out.dict
    )

    // 10. Aggregate Report
    // Mix all QC outputs into one channel for MultiQC
    FASTQC.out.zip
        .mix(TRIM_PRIMERS.out.log)
        .mix(PICARD_METRICS.out)
        .collect()
        .set { multiqc_inputs }

    MULTIQC(multiqc_inputs)
}