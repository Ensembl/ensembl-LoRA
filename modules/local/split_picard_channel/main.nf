process SPLIT_PICARD_CHANNEL {
    label 'process_very_low'
    label "error_retry"
    input:
        tuple path(bam), path(reference), val(genome_id), val(chromosome_id), val(chromosome_count), val(bin_count)
    output:
        path("picard_sample_sheet.csv"), emit: picard_sample_sheet
    script:
        """
        source ~/venv/bin/activate
        process_picard_output.py ${genome_id} ${reference} ${chromosome_id} ${chromosome_count} ${bin_count} picard_sample_sheet.csv
        """
}
