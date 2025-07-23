process COMBINE_FULL_LENGTH_READS {
    label 'process_medium'
    label 'error_retry'
    container "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    publishDir "${params.outDir}/${genome_id}/full_length/reads/merged/", mode: 'copy'
    input:
        tuple val(genome_id), path(bam), path(contains)
    output:
        tuple path("${params.runID}_${genome_id}_full_length.bam"), path("${params.runID}_${genome_id}_full_length.bam.csi"), emit:merged_full_length_bam
    script:
    """
    samtools merge -@ 20 -o ${params.runID}_${genome_id}_full_length.bam ${bam}
    samtools index -c -@ 20 ${params.runID}_${genome_id}_full_length.bam
    """
}