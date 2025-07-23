process FASTQ_TO_FASTA {
    label "process_very_low"
    label "error_retry"
    container "docker://quay.io/biocontainers/seqtk:1.4--he4a0461_1"
    input:
        tuple path(fastq), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species)
    output:
        tuple path("*.fasta.gz"), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species), emit: converted_fasta
    script:
        """
        fastq_to_fasta.sh ${fastq} ${fastq.simpleName}.fasta
        """
}