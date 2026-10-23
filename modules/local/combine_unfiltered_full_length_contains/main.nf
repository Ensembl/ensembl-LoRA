process COMBINE_UNFILTERED_FULL_LENGTH_CONTAINS {
    label 'process_very_low'
    label 'error_retry'
    input:
        tuple val(genome_id), path(full_length_contains)
    output:
        path("${params.runID}_${genome_id}_unfiltered_contains.tsv"), emit:example_ch
    when:
        task.ext.when == null || task.ext.when

    script:
    """
    source ~/venv/bin/activate
    combine_tmerge.py -o ${params.runID}_${genome_id}
    """
}