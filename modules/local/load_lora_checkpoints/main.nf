process LOAD_LORA_CHECKPOINTS {
    label "process_very_low"
    label "error_retry"

    input:
        val checkpoint_outdir
        path sample_metadata
        path reference_sheet

    output:
        path("checkpoint_plan.csv"), emit: checkpoint_plan
        path("read_prep_fasta.csv"), emit: read_prep_fasta
        path("mapped_bam.csv"), emit: mapped_bam
        path("canonical_bam.csv"), emit: canonical_bam
        path("chromosome_models.csv"), emit: chromosome_models
        path("prefiltered_models.csv"), emit: prefiltered_models
        path("prefiltered_contains.csv"), emit: prefiltered_contains
        path("filtered_models.csv"), emit: filtered_models
        path("filtered_chromosome_models.csv"), emit: filtered_chromosome_models
        path("full_length_stats.csv"), emit: full_length_stats
        path("read_summaries.csv"), emit: read_summaries
        path("sqanti3_dashboard.csv"), emit: sqanti3_dashboard
        path("sample_counts.csv"), emit: sample_counts
        path("input_read_counts.csv"), emit: input_read_counts
        path("prefiltered_model_counts.csv"), emit: prefiltered_model_counts
        path("postfiltered_model_counts.csv"), emit: postfiltered_model_counts
        path("prefiltered_sqanti3_counts.csv"), emit: prefiltered_sqanti3_counts
        path("postfiltered_sqanti3_counts.csv"), emit: postfiltered_sqanti3_counts
        path("data_checks.csv"), emit: data_checks
        path("dashboard_bam_stats_inputs.csv"), emit: dashboard_bam_stats_inputs

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        source ~/venv/bin/activate
        lora_checkpoint_manifest.py \\
            --outdir ${checkpoint_outdir} \\
            --run_id ${params.runID} \\
            --resume_from ${params.resume_from} \\
            --resume_scope ${params.resume_scope} \\
            --sample_metadata ${sample_metadata} \\
            --reference_sheet ${reference_sheet} \\
            --output_dir .
        """
}
