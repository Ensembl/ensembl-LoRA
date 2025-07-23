process COMBINE_FULL_LENGTH_STATS {
    label 'process_very_low'
    label 'error_retry'
    publishDir "${params.outDir}/${genome_id}/full_length/stats/", mode: 'copy', pattern: "*_stats.tsv"
    input:
        tuple val(genome_id), path(stats)
    output:
        path("${params.runID}_${genome_id}_full_length_stats.tsv"), emit: full_length_stats
    script:
    """
    combine_full_length_stats.py -f *_stats.tsv -p ${params.runID}_${genome_id}
    """
} 