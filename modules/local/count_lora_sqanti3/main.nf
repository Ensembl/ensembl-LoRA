process COUNT_LORA_SQANTI3 {
    label "process_medium"
    label "error_retry"

    input:
        tuple val(stage), val(genome_id), val(sqanti3_available), path(sqanti3_dir)

    output:
        tuple val(stage), val(genome_id), path("${params.runID}_${genome_id}_${stage}_sqanti3_counts.tsv"), emit: sqanti3_counts

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        source ~/venv/bin/activate
        count_lora_dashboard_metrics.py sqanti3 \\
            --run_id ${params.runID} \\
            --genome_id ${genome_id} \\
            --stage ${stage} \\
            --sqanti3_dir ${sqanti3_dir} \\
            --run_sqanti3 ${sqanti3_available} \\
            --output_dir .
        """
}
