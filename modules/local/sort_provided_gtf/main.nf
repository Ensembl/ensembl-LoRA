process SORT_PROVIDED_GTF {
    label "process_very_high"
    label "error_retry"
    input:
    tuple path(provided_gtf), val(genome_id), val(individual_bin_count), val(total_bin_count)
    output:
    tuple path ("*.sorted.gtf"), val(genome_id), val(total_bin_count), emit: sorted_gtf
    script:
    """
    sort -V -k1,1 -k4,4n ${provided_gtf} > ${provided_gtf.simpleName}_gtf.sorted.gtf
    """
}