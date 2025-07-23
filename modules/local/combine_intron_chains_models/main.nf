process COMBINE_INTRON_CHAINS_MODELS {
    label 'process_very_low'
    label 'error_retry'
    input:
        tuple val(genome_id), path(intron_chains), path(single_exons)
    output:
        tuple val(genome_id), path("${genome_id}.models_intron_chains.tsv"), path("${genome_id}.models_single_exons.tsv"), val(2), emit: merged_intron_chains_models
    script:
        """
        cat ${intron_chains} > ${genome_id}.models_intron_chains.tsv
        cat ${single_exons} > ${genome_id}.models_single_exons.tsv
        """
}
