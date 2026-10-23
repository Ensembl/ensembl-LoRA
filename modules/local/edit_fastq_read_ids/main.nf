process EDIT_FASTQ_IDS {
    label "process_low"
    label "multi_cpu_low"
    label "error_retry"
    input:
        tuple path(fastq), val(format), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species), val(read_id_map), val(original_reads), val(tissue), val(cell_line)
    output:
        tuple path("*.reads.fasta.gz"), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species), emit: edited_fasta
        tuple path("*.read_id_map.tsv.gz"), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species), emit: read_id_map
        tuple path("*.read_summary.tsv"), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species), emit: read_summary
    when:
        task.ext.when == null || task.ext.when

    script:
        """
        edit_fastq_ids.sh "${fastq}" "${read_id_map}" "${original_reads}" "${seq_source}" "${genome_id}" "${tissue}" "${cell_line}" "${task.cpus}" "${params.read_prep_compression_level}"
        """
}
