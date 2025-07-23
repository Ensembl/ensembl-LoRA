process SPLIT_GTF_INTO_BINS {
    label "process_medium"
    label "error_retry"
    input:
    tuple path(gtf), val(genome_id), val(individual_bin_count), val(total_bin_count)
    output:
    tuple path("${params.runID}_*.fixed.gtf"), val(genome_id), val(individual_bin_count), val(total_bin_count), emit: binned_gtf
    script:
    """
    split_gtf_into_bins.py ${gtf} ${individual_bin_count} ${params.runID}
    """
}