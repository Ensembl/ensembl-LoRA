process BUILD_LORA_DATA_CHECKS {
    label "process_medium"
    label "error_retry"

    input:
        path sample_metadata
        val sample_counts
        val input_read_counts
        val model_counts
        val sqanti3_counts
        val full_length_stats
        val bam_stats

    output:
        path("data_checks/*"), emit: data_checks

    when:
        task.ext.when == null || task.ext.when

    script:
        def sampleCountsArgs = shellArgs(sample_counts)
        def inputReadCountsArgs = shellArgs(input_read_counts)
        def modelCountsArgs = shellArgs(model_counts)
        def sqanti3CountsArgs = shellArgs(sqanti3_counts)
        def fullLengthStatsArgs = shellArgs(full_length_stats)
        def bamStatsArgs = shellArgs(bam_stats)
        """
        source ~/venv/bin/activate
        build_lora_data_checks.py \\
            --run_id ${params.runID} \\
            --sample_metadata ${sample_metadata} \\
            --sample_counts ${sampleCountsArgs} \\
            --input_read_counts ${inputReadCountsArgs} \\
            --model_counts ${modelCountsArgs} \\
            --sqanti3_counts ${sqanti3CountsArgs} \\
            --full_length_stats ${fullLengthStatsArgs} \\
            --bam_stats ${bamStatsArgs} \\
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
