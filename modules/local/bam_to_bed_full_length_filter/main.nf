process BAM_TO_BED_FULL_LENGTH_FILTER {
    label "process_very_low"
    label "error_retry"
    container "docker://quay.io/biocontainers/bedtools:2.24--1"
    input:
        tuple path(bam), path(bai), val(genome_id), val(chromosome_id), val(chromosome_count)
    output:
        tuple path("*.bed"), val(genome_id), val(chromosome_id), val(chromosome_count), emit: converted_bed
    script:
        """
        bedtools bamtobed -bed12 -i ${bam} > ${bam.simpleName}.bed
        """
}