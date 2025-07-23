process PROCESS_SAMPLE_SHEETS {
    label "process_very_low"
    label "error_retry"
    input:
        path read_sample_sheet
        path reference_sample_sheet
    output:
        path("processed_sample_sheets.csv") , emit: processed_sample_sheets
    script:
        """
        source ~/venv/bin/activate
        process_sample_sheets.py ${read_sample_sheet} ${reference_sample_sheet} processed_sample_sheets.csv
        """
}