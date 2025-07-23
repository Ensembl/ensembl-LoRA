process INTRON_CHAIN_EXTRACTOR_READS {
    label 'process_low'
    label 'error_retry'
    input:
        tuple path(bed), val(genome_id), val(chromosome_id), val(chromosome_count)
    output:
        tuple val(chromosome_id), val(genome_id), val(chromosome_count), path("${params.runID}_${chromosome_id}_reads_intron_chains.tsv"), path("${params.runID}_${chromosome_id}_reads_single_exons.tsv"), emit: intron_chains_reads
    script:
        """
        source ~/venv/bin/activate
        extract_intron_chains_from_bed12.py ${bed} ${params.runID}_${chromosome_id}_reads
        """
}
