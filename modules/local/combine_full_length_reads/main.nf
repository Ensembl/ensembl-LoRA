process COMBINE_FULL_LENGTH_READS {
    label 'process_medium'
    label 'multi_cpu_high'
    label 'error_retry'
    container "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    input:
        tuple val(genome_id), path(bam), path(contains)
    output:
        tuple path("${params.runID}_${genome_id}_full_length.bam"), path("${params.runID}_${genome_id}_full_length.bam.csi"), emit:merged_full_length_bam
    when:
        task.ext.when == null || task.ext.when

    script:
    """
    samtools merge -@ ${task.cpus} -o ${params.runID}_${genome_id}_full_length.bam ${bam}
    samtools index -c -@ ${task.cpus} ${params.runID}_${genome_id}_full_length.bam
    """
}
