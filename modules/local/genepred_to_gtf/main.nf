process GENEPRED_TO_GTF {
    label "process_very_low"
    label "error_retry"
    container "docker://quay.io/biocontainers/ucsc-genepredtogtf:377--ha8a8165_5"
    input:
        tuple path(genePred), val(genome_id), val(chromosome_id), val(chromosome_count), val(bin_count)
    output:
        tuple path("*.gtf"), val(genome_id), val(chromosome_id), val(chromosome_count), val(bin_count), emit: converted_gtf
    script:
        """
        genePredToGtf file ${genePred} ${genePred.simpleName}.gtf
        """
}