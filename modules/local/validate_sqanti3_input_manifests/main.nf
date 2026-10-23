process VALIDATE_SQANTI3_INPUT_MANIFESTS {
    label "process_very_low"
    label "error_retry"

    input:
        path(input_manifests, stageAs: "sqanti3_input_manifests/*")

    output:
        path("sqanti3_qc_inputs.tsv"), emit: sqanti3_qc_inputs

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        source ~/venv/bin/activate
        validate_sqanti3_input_manifests.py \\
            --manifest_root sqanti3_input_manifests \\
            --output sqanti3_qc_inputs.tsv
        """
}
