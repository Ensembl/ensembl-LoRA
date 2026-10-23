process COMBINE_FULL_LENGTH_STATS {
    label 'process_very_low'
    label 'error_retry'
    input:
        tuple val(genome_id), path(stats)
    output:
        tuple val(genome_id), path("${params.runID}_${genome_id}_full_length_stats.tsv"), emit: full_length_stats
    when:
        task.ext.when == null || task.ext.when

    script:
    """
    combine_full_length_stats.py -f *_stats.tsv -p ${params.runID}_${genome_id}
    """
}
