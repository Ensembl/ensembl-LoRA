process MERGE_SQANTI3_OUTPUTS {
    label "process_medium"
    label "error_retry"

    input:
        tuple val(stage), val(genome_id), val(chunk_ids), path(chunk_manifests, stageAs: "sqanti3_chunk_manifests/*")

    output:
        tuple val(stage), val(genome_id), path("sqanti3_output/*"), emit: sqanti3_results

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        source ~/venv/bin/activate
        find sqanti3_chunk_manifests -maxdepth 2 -type f -name '*.sqanti3_chunk_manifest.tsv' | sort

        merge_sqanti3_outputs.py \\
            --manifest_root sqanti3_chunk_manifests \\
            --output_dir sqanti3_output \\
            --run_id ${params.runID} \\
            --genome_id ${genome_id} \\
            --stage ${stage}
        """
}
