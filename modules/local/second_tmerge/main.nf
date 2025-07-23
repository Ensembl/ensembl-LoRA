process SECOND_TMERGE {
    label 'process_low'
    label 'error_retry'
    input:
        tuple val(genome_id), path(merged_gtf), path(initial_merged_contains), val(chromosome_count), val(chromosome_id)
    output:
        tuple val(genome_id), path("${merged_gtf.simpleName}.unfixed.gtf"), path(initial_merged_contains), val(chromosome_count), val(chromosome_id), emit: second_tmerge_output
    script:
        """
        tmerge.pl --exonOverhangTolerance 10 ${merged_gtf} > ${merged_gtf.simpleName}.unfixed.gtf
        """
}
