process COUNT_LORA_MODELS {
    label "process_medium"
    label "error_retry"

    input:
        tuple val(stage), val(genome_id), path(gtf), path(contains)

    output:
        tuple val(stage), val(genome_id), path("${params.runID}_${genome_id}_${stage}_model_counts.tsv"), emit: model_counts

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        source ~/venv/bin/activate
        count_lora_dashboard_metrics.py models \\
            --run_id ${params.runID} \\
            --genome_id ${genome_id} \\
            --stage ${stage} \\
            --gtf ${gtf} \\
            --contains ${contains} \\
            --output_dir .
        """
}
