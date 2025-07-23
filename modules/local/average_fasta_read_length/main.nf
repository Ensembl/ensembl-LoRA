process AVERAGE_FASTA_READ_LENGTH {
    label "process_very_low"
    label "error_retry"
    input:
        tuple path(fasta), val(format), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species)
    output:
        tuple path("*_fasta.gz"), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species), optional: true, emit: passed
        path("*.failed"), optional: true , emit: failed
    script:
        """
        average_fasta_read_length.sh ${fasta}
        """
}
