process SUMMARISE_LORA_DASHBOARD {
    label "process_medium"
    label "error_retry"

    input:
        tuple val(genome_id), path(dashboard_inputs)

    output:
        tuple val(genome_id),
            path("${params.runID}_${genome_id}_run_summary.json"),
            path("${params.runID}_${genome_id}_dashboard_data.json"),
            emit: dashboard_summary

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        source ~/venv/bin/activate
        summarise_lora_dashboard.py \\
            --run_id ${params.runID} \\
            --genome_id ${genome_id} \\
            --dashboard_input_dir ${dashboard_inputs} \\
            --output_dir . \\
            --min_fl ${params.minFL} \\
            --genome_size ${params.genomeSize}
        """
}
