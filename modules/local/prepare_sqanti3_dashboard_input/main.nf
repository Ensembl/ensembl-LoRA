process PREPARE_SQANTI3_DASHBOARD_INPUT {
    label "process_medium"
    label "error_retry"

    input:
        tuple val(genome_id), path(prefiltered_sqanti3_files, stageAs: "prefiltered_sqanti3/*"), path(postfiltered_sqanti3_files, stageAs: "postfiltered_sqanti3/*")

    output:
        tuple val(genome_id), val(true), path("sqanti3"), emit: sqanti3_dashboard_input

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        mkdir -p sqanti3/prefiltered sqanti3/postfiltered
        cp -R prefiltered_sqanti3/. sqanti3/prefiltered/
        cp -R postfiltered_sqanti3/. sqanti3/postfiltered/
        """
}
