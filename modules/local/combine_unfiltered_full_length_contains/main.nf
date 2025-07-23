process COMBINE_UNFILTERED_FULL_LENGTH_CONTAINS {
    label 'process_very_low'
    label 'error_retry'
    publishDir "${params.outDir}/${genome_id}/full_length/all_full_length/", mode: 'copy'
    input:
        tuple val(genome_id), path(full_length_contains)
    output:
        path("${params.runID}_${genome_id}_unfiltered_contains.tsv"), emit:example_ch
    script:
    """
    source ~/venv/bin/activate
    combine_tmerge.py -o ${params.runID}_${genome_id}
    """
}