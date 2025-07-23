process COMBINE_FILTERED_MODELS {
    label 'process_very_low'
    label 'error_retry'
    publishDir "${params.outDir}/${genome_id}/full_length/models/merged/", mode: 'copy'
    input:
        tuple val(genome_id), path(gtf), path(full_length_contains)
    output:
        tuple path("*_filtered_transcripts.gtf"), path("*_filtered_contains.tsv"), emit:filtered_gtf
    script:
    """
    source ~/venv/bin/activate
    combine_tmerge.py -r filter -o ${params.runID}_${genome_id}
    """
}