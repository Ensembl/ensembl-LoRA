process AVERAGE_FASTQ_READ_LENGTH {
    label "process_very_low"
    label "error_retry"
    input:
        tuple path(fastq), val(format), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species), val(read_id_map), val(original_reads), val(tissue), val(cell_line)
    output:
        tuple path("*_fastq.gz"), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species), val(read_id_map), val(original_reads), val(tissue), val(cell_line), optional: true, emit: passed
        path("*.failed"), optional: true , emit: failed
    when:
        task.ext.when == null || task.ext.when

    script:
        """
        average_fastq_read_length.sh ${fastq}
        """
}
