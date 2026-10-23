process COMBINE_FILTERED_MODELS {
    label 'process_very_low'
    label 'error_retry'
    input:
        tuple val(genome_id), path(gtf), path(full_length_contains)
    output:
        tuple val(genome_id), path("*_filtered_transcripts.gtf"), path("*_filtered_contains.tsv"), emit:filtered_gtf
    when:
        task.ext.when == null || task.ext.when

    script:
    """
    source ~/venv/bin/activate
    combine_tmerge.py -r filter -o ${params.runID}_${genome_id}
    """
}
