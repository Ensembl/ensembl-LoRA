process FULL_LENGTH_CONTAINS {
    label 'process_low'
    label 'error_retry'
    input:
        tuple val(chromosome_id), val(genome_id), val(chromosome_count), path(intron_chains_reads), path(single_exons_reads), path(intron_chains_models), path(single_exons_models), path(contains_file)
    output:
        tuple val(chromosome_id), val(genome_id), path("${params.runID}_${chromosome_id}_full_length_contains.tsv"), path("${params.runID}_${chromosome_id}_stats.tsv"), val(chromosome_count), emit: full_length_contains
    script:
        """
        source ~/venv/bin/activate
        full_length_read_extraction.py ${intron_chains_models} ${single_exons_models} ${intron_chains_reads} ${single_exons_reads} ${contains_file} ${params.runID}_${chromosome_id}
        """
}
