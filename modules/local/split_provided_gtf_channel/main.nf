process SPLIT_PROVIDED_GTF_CHANNEL {
    label 'process_very_low'
    label 'error_retry'
    input:
        tuple path(binned_gtfs), val(genome_id), val(individual_bin_count), val(total_bin_count)
    output:
        path("prepared_gtf_sample_sheet.csv"), emit: prepared_gtf_sample_sheet
    script:
        """
        source ~/venv/bin/activate
        process_provided_gtf_split.py ${genome_id} ${individual_bin_count} ${total_bin_count} prepared_gtf_sample_sheet.csv
        """
}
