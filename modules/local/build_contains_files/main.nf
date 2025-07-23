process BUILD_CONTAINS_FILES {
    label "process_very_low"
    label "error_retry"
    publishDir "${params.outDir}/${genome_id}/all_models/chromosomes/", mode: 'copy'
    input:
        tuple val(genome_id), path(second_tmerge_gtf), path(initial_merged_contains), path(final_merged_contains), val(chromosome_count), val(chromosome_id)
    output:
        tuple val(genome_id), path("${params.runID}_${chromosome_id}.combined.gtf"), path("${params.runID}_${chromosome_id}.combined.contains.tsv"), val(chromosome_count), val(chromosome_id), emit:combined_contains
    script:
        """
        if [ -s ${initial_merged_contains} ]; then
            source ~/venv/bin/activate
            combine_tmerge_contains_outputs.py ${initial_merged_contains} ${final_merged_contains} ${params.runID}_${chromosome_id}.combined
            cp ${second_tmerge_gtf} ${params.runID}_${chromosome_id}.combined.gtf
        else
            touch ${params.runID}_${chromosome_id}.combined.gtf
            touch ${params.runID}_${chromosome_id}.combined.contains.tsv
        fi
        """
}
