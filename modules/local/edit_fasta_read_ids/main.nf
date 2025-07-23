process EDIT_FASTA_IDS {
    label "process_very_low"
    label "error_retry"
    input:
        tuple path(fasta), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species)
    output:
        tuple path("*.reads.fasta.gz"), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species), emit: edited_fasta
    script:
        """
        edit_fasta_ids.sh ${fasta}
        """
}
