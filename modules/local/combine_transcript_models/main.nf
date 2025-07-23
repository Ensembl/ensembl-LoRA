process COMBINE_TRANSCRIPT_MODELS {
    label "process_very_low"
    label "error_retry"
    publishDir "${params.outDir}/${genome_id}/all_models/merged/", mode: 'copy'
    input:
        tuple val(genome_id), path(final_chromosome_gtf), path(final_chromosome_contains), val(chromosome_id)
    output:
        tuple val(genome_id), path("${params.runID}_${genome_id}.gtf"), emit: combined_transcript_models
        tuple val(genome_id), path("${params.runID}_${genome_id}_contains.tsv"), emit: combined_contains_files
    script:
        """
        source ~/venv/bin/activate
        combine_tmerge.py -r final -o ${params.runID}_${genome_id}
        """
}
