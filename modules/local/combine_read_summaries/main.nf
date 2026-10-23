process COMBINE_READ_SUMMARIES {
    label 'process_medium'
    label 'error_retry'

    input:
        tuple val(genome_id), path(read_summaries)

    output:
        tuple val(genome_id), path("${params.runID}_${genome_id}_read_summary.tsv"), emit: combined_read_summary

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    awk 'FNR == 1 && NR != 1 { next } { print }' *.read_summary.tsv > ${params.runID}_${genome_id}_read_summary.tsv
    """
}
