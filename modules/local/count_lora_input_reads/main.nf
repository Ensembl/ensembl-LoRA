process COUNT_LORA_INPUT_READS {
    label "process_medium"
    label "error_retry"

    input:
        tuple val(genome_id), path(read_summary)

    output:
        tuple val(genome_id), path("${params.runID}_${genome_id}_input_read_counts.tsv"), emit: input_read_counts

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        source ~/venv/bin/activate
        count_lora_dashboard_metrics.py input_reads \\
            --run_id ${params.runID} \\
            --genome_id ${genome_id} \\
            --read_summary ${read_summary} \\
            --output_dir .
        """
}
