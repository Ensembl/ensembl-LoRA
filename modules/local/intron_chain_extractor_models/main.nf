process INTRON_CHAIN_EXTRACTOR_MODELS {
    label 'process_low'
    label 'error_retry'
    input:
        tuple val(genome_id), path(gtf), path(contains), val(chromosome_count), val(chromosome_id)
    output:
        tuple val(chromosome_id), val(genome_id), val(chromosome_count), path("${params.runID}_${chromosome_id}_models_intron_chains.tsv"), path("${params.runID}_${chromosome_id}_models_single_exons.tsv"), path(contains), emit: intron_chains_models
    script:
        """
        if [ -s ${gtf} ]; then
            source ~/venv/bin/activate
            proto_transcript_introns.py -t ${gtf} -p ${params.runID}_${chromosome_id}_models
        else
            touch ${params.runID}_${chromosome_id}_models_intron_chains.tsv
            touch ${params.runID}_${chromosome_id}_models_single_exons.tsv
        fi
        """
}
