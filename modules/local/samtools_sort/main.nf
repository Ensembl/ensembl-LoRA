process SAMTOOLS_SORT {
    label "process_high"
    label "multi_cpu_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    input:
        tuple path(bam), path(reference_genome), val(genome_id), val(raw_file_count)
    output: 
        tuple val(genome_id), path("*.bam"), path(reference_genome), val(raw_file_count), emit: sorted_bam
        tuple val(genome_id), path("*.bam.csi"), path(reference_genome), val(raw_file_count), emit: sorted_bam_bai
        tuple val(genome_id), path("*.bam"), path("*.bam.csi"), path(reference_genome), val(raw_file_count), emit: sorted_bam_with_index
    when:
        task.ext.when == null || task.ext.when

    script:
        """
        samtools sort -@ ${task.cpus} -O bam -o ${bam.simpleName}_${genome_id}.sorted.bam ${bam}
        samtools index -c -@ ${task.cpus} ${bam.simpleName}_${genome_id}.sorted.bam
        if [ \$? -eq 0 ]; then
            rm -f ${bam} \$(readlink -f ${bam})
        fi
        """
}
