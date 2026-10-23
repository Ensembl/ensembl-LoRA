process SAMTOOLS_MERGE {
    label "process_medium"
    label "multi_cpu_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    input:
        tuple path(bams), path(reference_genome), val(genome_id)
    output:
        tuple path("${params.runID}_${genome_id}.bam"), path("${params.runID}_${genome_id}.bam.csi"), path(reference_genome), val(genome_id), emit: merged_bam
    when:
        task.ext.when == null || task.ext.when

    script:
        """
        samtools merge -@ ${task.cpus} -o ${params.runID}_${genome_id}.bam ${bams}
        samtools index -c -@ ${task.cpus} ${params.runID}_${genome_id}.bam
        if [ \$? -eq 0 ]; then
            rm -f ${bams} \$(readlink -f ${bams})
        fi
        """
}
