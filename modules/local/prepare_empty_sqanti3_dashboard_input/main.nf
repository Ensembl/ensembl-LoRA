process PREPARE_EMPTY_SQANTI3_DASHBOARD_INPUT {
    label "process_medium"
    label "error_retry"

    input:
        val genome_id

    output:
        tuple val(genome_id), val(false), path("sqanti3"), emit: sqanti3_dashboard_input

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        mkdir -p sqanti3/prefiltered sqanti3/postfiltered
        """
}
