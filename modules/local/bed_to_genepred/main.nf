process BED_TO_GENEPRED {
    label "process_very_low"
    label "error_retry"
    container "docker://quay.io/biocontainers/ucsc-bedtogenepred:377--ha8a8165_3"
    input:
        tuple  path(bed), val(genome_id), val(chromosome_id), val(chromosome_count), val(bin_count)
    output:
        tuple path("*.txt"), val(genome_id), val(chromosome_id), val(chromosome_count), val(bin_count), emit: converted_genepred
    script:
        """
        bedToGenePred ${bed} ${bed.simpleName}.txt
        """
}