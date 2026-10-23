process WRITE_SQANTI3_INPUT_MANIFEST {
    label "process_very_low"
    label "error_retry"

    input:
        tuple val(stage),
            val(genome_id),
            val(chunk_id),
            path(isoforms),
            path(reference_fasta),
            path(reference_gtf),
            path(polya_peak),
            val(polya_source),
            path(cage_peak),
            val(cage_source)

    output:
        path("*.sqanti3_input_manifest.tsv"), emit: sqanti3_input_manifest

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = "${params.runID}_${genome_id}_${stage}_${chunk_id}"
        """
        manifest="\$(mktemp ${prefix}.XXXXXX.sqanti3_input_manifest.tsv)"
        {
            printf 'stage\\tgenome_id\\tchunk_id\\tisoforms\\treference_fasta\\treference_gtf\\tpolya_peak\\tpolya_source\\tcage_peak\\tcage_source\\n'
            printf '${stage}\\t${genome_id}\\t${chunk_id}\\t%s\\t%s\\t%s\\t%s\\t${polya_source}\\t%s\\t${cage_source}\\n' \\
                "\$(realpath ${isoforms})" \\
                "\$(realpath ${reference_fasta})" \\
                "\$(realpath ${reference_gtf})" \\
                "\$(realpath ${polya_peak})" \\
                "\$(realpath ${cage_peak})"
        } >| "\$manifest"
        """
}
