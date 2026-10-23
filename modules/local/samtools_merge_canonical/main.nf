process SAMTOOLS_MERGE_CANONICAL {
    label "process_medium"
    label "multi_cpu_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    input:
        tuple val(genome_id), path(bams), path(bai)
    output:
        tuple path("${params.runID}_${genome_id}.canonical.bam"), path("${params.runID}_${genome_id}.canonical.bam.csi"), val(genome_id), emit: merged_bam
    when:
        task.ext.when == null || task.ext.when

    script:
        """
        samtools merge -@ ${task.cpus} -o ${params.runID}_${genome_id}.canonical.bam ${bams}
        samtools index -c -@ ${task.cpus} ${params.runID}_${genome_id}.canonical.bam
        """
}
