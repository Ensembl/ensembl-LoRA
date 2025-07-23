process COMBINE_FIRST_TMERGE {
    label "process_very_low"
    label "error_retry"
    input:
    tuple val(genome_id), path(gtfs), path(contains), val(chromosome_count), val(chromosome_id)
    output:
    tuple val(genome_id), path("${params.runID}_${chromosome_id}_initial_transcripts.gtf"), path("${params.runID}_${chromosome_id}_initial_contains.tsv"), val(chromosome_count), val(chromosome_id), emit: combined_bins
    script:
    """
    source ~/venv/bin/activate
    combine_tmerge.py -r first -o ${params.runID}_${chromosome_id}
    """
}
