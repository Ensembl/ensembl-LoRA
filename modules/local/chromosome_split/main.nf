process CHROMOSOME_SPLIT {
    label "process_very_low"
    label "multi_cpu_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    input:
        tuple path(bam), path(bai), path(reference_genome), val(genome_id), val(chromosome), val(chromosome_id), val(chromosome_count), val(bin_count) 
    output:
        tuple path("${bam.simpleName}_${chromosome}.bam"), path("${bam.simpleName}_${chromosome}.bam.csi"), path(reference_genome), val(genome_id), val(chromosome_id), val(chromosome_count), val(bin_count), emit: chromosome_bam
    when:
        task.ext.when == null || task.ext.when

    script:
        """
        samtools view -@ ${task.cpus} -b ${bam} ${chromosome} | samtools sort -@ ${task.cpus} -O bam -o ${bam.simpleName}_${chromosome}.bam
        samtools index -c -@ ${task.cpus} ${bam.simpleName}_${chromosome}.bam
        """
}
