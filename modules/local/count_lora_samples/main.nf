process COUNT_LORA_SAMPLES {
    label "process_medium"
    label "error_retry"

    input:
        tuple val(genome_id), path(sample_metadata)

    output:
        tuple val(genome_id), path("${params.runID}_${genome_id}_sample_counts.tsv"), emit: sample_counts

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        source ~/venv/bin/activate
        count_lora_dashboard_metrics.py samples \\
            --run_id ${params.runID} \\
            --genome_id ${genome_id} \\
            --sample_metadata ${sample_metadata} \\
            --output_dir .
        """
}
