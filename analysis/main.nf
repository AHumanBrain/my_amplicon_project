#!/usr/bin/env nextflow

// --- Input Parameters ---
params.reads = "${projectDir}/../raw_data/*_R{1,2}.fastq.gz"
params.genome = "${projectDir}/../common_refs/ecoli_genome.fna"
params.primers = "${projectDir}/../design/output/panel_pool_1_trimming.fasta"
params.bed = "${projectDir}/../design/output/panel_pool_1.bed"
params.outdir = "${projectDir}/../results"
// HaplotypeCaller assembly settings tuned for amplicon data (see README section 4):
// every read starts at the primer, so SNVs near read ends have no left/right
// anchor and clustered SNVs defeat the default k-mers (10,25).
params.hc_args = "-kmer-size 10 -kmer-size 15 -kmer-size 25 --recover-all-dangling-branches --min-dangling-branch-length 1"

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

    stub:
    """
    touch fastqc_stub.html fastqc_stub.zip
    """

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

    stub:
    """
    touch ${sample_id}_R1.trimmed.fastq.gz ${sample_id}_R2.trimmed.fastq.gz ${sample_id}_cutadapt.log
    """

    // The trimming FASTA holds <target>_F and <target>_R oligos. Reads begin at the
    // primer, so each primer is matched anchored at the 5' end (^file:), and R1 is only
    // searched for forward primers and R2 for reverse primers.
    script:
    def common = "--discard-untrimmed -j ${task.cpus}"
    if (reads instanceof List && reads.size() == 2) {
        """
        awk '/^>/{p=(\$0 ~ /_F\$/)} p' ${primers} > fwd_primers.fa
        awk '/^>/{p=(\$0 ~ /_R\$/)} p' ${primers} > rev_primers.fa
        cutadapt ${common} \
            -g ^file:fwd_primers.fa \
            -G ^file:rev_primers.fa \
            -o ${sample_id}_R1.trimmed.fastq.gz \
            -p ${sample_id}_R2.trimmed.fastq.gz \
            ${reads[0]} ${reads[1]} > ${sample_id}_cutadapt.log
        """
    } else {
        """
        awk '/^>/{p=(\$0 ~ /_F\$/)} p' ${primers} > fwd_primers.fa
        cutadapt ${common} \
            -g ^file:fwd_primers.fa \
            -o ${sample_id}.trimmed.fastq.gz \
            ${reads} > ${sample_id}_cutadapt.log
        """
    }
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

    stub:
    """
    touch ${genome}.amb ${genome}.ann ${genome}.bwt ${genome}.pac ${genome}.sa
    touch ${genome}.fai
    touch ${genome.baseName}.dict
    """

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
    tuple val(sample_id), path("*.sorted.bam"), emit: bam

    stub:
    """
    touch ${sample_id}.sorted.bam
    """

    script:
    if (reads instanceof List && reads.size() == 2) {
        """
        bwa mem -t ${task.cpus} -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA' $genome ${reads[0]} ${reads[1]} | \
        samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
        """
    } else {
        """
        bwa mem -t ${task.cpus} -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA' $genome ${reads} | \
        samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
        """
    }
}

process SORT_INDEX_BAM {
    // BAMs arrive coordinate-sorted from ALIGN; this step only indexes them.
    tag "$sample_id"
    publishDir "${params.outdir}/bams", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), path("*.sorted.bam.bai"), emit: sorted_bam

    stub:
    """
    touch ${sample_id}.sorted.bam.bai
    """

    script:
    """
    samtools index $bam
    """
}

process PREP_INTERVALS {
    // Converts simple BED file to Picard IntervalList
    input:
    path bed
    path dict

    output:
    path "targets.interval_list", emit: interval_list

    stub:
    """
    touch targets.interval_list
    """

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

    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path fai
    path dict
    path interval_list

    output:
    path "${sample_id}_pcr_metrics.txt", emit: pcr_metrics
    path "${sample_id}.per_target_coverage.txt", emit: target_coverage
    path "*.{txt,pdf}", emit: all_metrics

    stub:
    """
    touch ${sample_id}_pcr_metrics.txt ${sample_id}.per_target_coverage.txt \
          ${sample_id}_alignment_metrics.txt ${sample_id}_insert_size_metrics.txt
    """

    script:
    """
    # 1. Alignment Metrics
    picard CollectAlignmentSummaryMetrics \
    R=$fasta I=$bam O=${sample_id}_alignment_metrics.txt

    # 2. Insert Size Metrics
    # Insert size needs paired reads; tolerate its absence (e.g. single-end data)
    picard CollectInsertSizeMetrics \
    I=$bam O=${sample_id}_insert_size_metrics.txt H=${sample_id}_insert_size_histogram.pdf || true

    picard CollectTargetedPcrMetrics \
    I=$bam O=${sample_id}_pcr_metrics.txt \
    AMPLICON_INTERVALS=$interval_list \
    TARGET_INTERVALS=$interval_list \
    PER_TARGET_COVERAGE=${sample_id}.per_target_coverage.txt \
    R=$fasta
    """
}

process GERMLINE_CALLING {
    tag "$sample_id"
    publishDir "${params.outdir}/variants/germline", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path fai
    path dict

    output:
    path "${sample_id}.g.vcf.gz", emit: vcf

    stub:
    """
    touch ${sample_id}.g.vcf.gz
    """

    // Using PCR indel model conservatively for germline
    script:
    """
    gatk HaplotypeCaller \
    -R $fasta \
    -I $bam \
    -O ${sample_id}.g.vcf.gz \
    --pcr-indel-model CONSERVATIVE \
    ${params.hc_args} \
    -ERC GVCF
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path '*' // Collects all files from previous processes

    output:
    path "multiqc_report.html"

    stub:
    """
    touch multiqc_report.html
    """

    script:
    """
    multiqc .
    """
}

process REPORT {
    tag "$sample_id"
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path pcr_metrics
    path target_coverage

    output:
    path "${sample_id}_balancing_feedback.json", emit: json
    path "${sample_id}_target_summary.tsv", emit: tsv
    path "${sample_id}_uniformity.png", emit: png

    stub:
    """
    touch ${sample_id}_balancing_feedback.json ${sample_id}_target_summary.tsv ${sample_id}_uniformity.png
    """

    script:
    """
    python ${projectDir}/balance.py -m $pcr_metrics -c $target_coverage \\
        -j ${sample_id}_balancing_feedback.json -t ${sample_id}_target_summary.tsv
    python ${projectDir}/plot_uniformity.py -m $pcr_metrics -c $target_coverage \\
        -o ${sample_id}_uniformity.png -s "${sample_id}"
    """
}

// --- Workflow ---

workflow {
    read_pairs_ch = Channel.fromFilePairs(params.reads, size: -1) { file -> file.name.replaceAll(/_R[12].*/, '') }

    // Reference inputs must be file objects; .first() makes value channels reused by every sample.
    genome_ch  = Channel.fromPath(params.genome,  checkIfExists: true).first()
    bed_ch     = Channel.fromPath(params.bed,     checkIfExists: true).first()
    primers_ch = Channel.fromPath(params.primers, checkIfExists: true).first()

    BUILD_INDICES(genome_ch)
    PREP_INTERVALS(bed_ch, BUILD_INDICES.out.dict)
    FASTQC(read_pairs_ch)
    TRIM_PRIMERS(read_pairs_ch, primers_ch)
    ALIGN(TRIM_PRIMERS.out.trimmed_reads, BUILD_INDICES.out.bwa_index, genome_ch)
    SORT_INDEX_BAM(ALIGN.out.bam)

    PICARD_METRICS(
        SORT_INDEX_BAM.out.sorted_bam,
        BUILD_INDICES.out.fasta, BUILD_INDICES.out.fai, BUILD_INDICES.out.dict,
        PREP_INTERVALS.out.interval_list
    )

    // Variants are called on the sorted BAM, NOT a duplicate-marked one: amplicon reads share
    // start positions, so GATK's default duplicate filter would discard nearly everything.
    GERMLINE_CALLING(
        SORT_INDEX_BAM.out.sorted_bam,
        BUILD_INDICES.out.fasta, BUILD_INDICES.out.fai, BUILD_INDICES.out.dict
    )

    REPORT(SORT_INDEX_BAM.out.sorted_bam, PICARD_METRICS.out.pcr_metrics, PICARD_METRICS.out.target_coverage)

    FASTQC.out.zip
        .mix(TRIM_PRIMERS.out.log, PICARD_METRICS.out.all_metrics, REPORT.out.png)
        .collect()
        .set { multiqc_inputs }
    MULTIQC(multiqc_inputs)
}
