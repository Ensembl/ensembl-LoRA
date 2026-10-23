process BUILD_LORA_DASHBOARD {
    label "process_medium"
    label "error_retry"

    input:
        tuple val(genome_id), path(dashboard_data)

    output:
        tuple val(genome_id), path("index.html"), path("dashboard_data.json"), emit: dashboard

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        source ~/venv/bin/activate
        build_lora_dashboard.py \\
            --dashboard_data ${dashboard_data} \\
            --output_dir .
        """
}
