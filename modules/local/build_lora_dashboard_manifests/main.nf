process BUILD_LORA_DASHBOARD_MANIFESTS {
    label "process_medium"
    label "error_retry"

    input:
        val sample_counts
        val input_read_counts
        val model_counts
        val sqanti3_counts
        val full_length_stats
        val data_checks
        path sample_metadata
        path reference_sheet

    output:
        path("dashboard_inputs/*"), emit: dashboard_inputs

    when:
        task.ext.when == null || task.ext.when

    script:
        def sampleCountsArgs = shellArgs(sample_counts)
        def inputReadCountsArgs = shellArgs(input_read_counts)
        def modelCountsArgs = shellArgs(model_counts)
        def sqanti3CountsArgs = shellArgs(sqanti3_counts)
        def fullLengthStatsArgs = shellArgs(full_length_stats)
        def dataChecksArgs = shellArgs(data_checks)
        """
        source ~/venv/bin/activate
        build_lora_dashboard_manifests.py \\
            --run_id ${params.runID} \\
            --sample_metadata ${sample_metadata} \\
            --reference_sheet ${reference_sheet} \\
            --sample_counts ${sampleCountsArgs} \\
            --input_read_counts ${inputReadCountsArgs} \\
            --model_counts ${modelCountsArgs} \\
            --sqanti3_counts ${sqanti3CountsArgs} \\
            --full_length_stats ${fullLengthStatsArgs} \\
            --data_checks ${dataChecksArgs} \\
            --output_dir .
        """
}

def shellArgs(values) {
    def items = values instanceof Collection ? values : [values]
    return items
        .findAll { it != null && it.toString() }
        .collect { "'" + it.toString().replace("'", "'\"'\"'") + "'" }
        .join(' ')
}
