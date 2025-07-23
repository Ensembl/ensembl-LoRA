process EDIT_FASTQ_IDS {
    label "process_very_low"
    label "error_retry"
    input:
        tuple path(fastq), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species)
    output:
        tuple path("*.reads.fastq.gz"), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species), emit: edited_fastq
    script:
        """
        edit_fastq_ids.sh ${fastq} 
        """
}
