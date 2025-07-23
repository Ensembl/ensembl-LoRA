process SECOND_TMERGE_TO_GTF {
    label "process_very_low"
    label "error_retry"
    input:
        tuple val(genome_id), path(retmerged_gtf), path(initial_merged_contains), val(chromosome_count), val(chromosome_id)
    output:
        tuple val(genome_id), path("*.final.tama.gtf"), path(initial_merged_contains), path("*.final.contains.tsv"), val(chromosome_count), val(chromosome_id), emit: retmerge_to_gtf_output
    script:
        """
        if [ -s ${retmerged_gtf} ]; then
            source ~/venv/bin/activate
            tmerge_to_gtf.py ${retmerged_gtf} ${params.runID} ${retmerged_gtf.simpleName}.final
        else
            touch ${retmerged_gtf.simpleName}.final.tama.gtf
            touch ${retmerged_gtf.simpleName}.final.contains.tsv
        fi
        """
}
