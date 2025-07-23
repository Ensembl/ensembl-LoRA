process AVERAGE_FASTQ_READ_LENGTH {
    label "process_very_low"
    label "error_retry"
    input:
        tuple path(fastq), val(format), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species)
    output:
        tuple path("*_fastq.gz"), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species), optional: true, emit: passed
        path("*.failed"), optional: true , emit: failed
    script:
        """
        average_fastq_read_length.sh ${fastq}
        """
}
