process SQANTI3_QC {
    label "process_high"
    label "error_retry"

    container "docker://anaconesalab/sqanti3:v6.0.1"

    cpus { params.sqanti3_cpus as int }

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
        tuple val(stage), val(genome_id), val(chunk_id), path("*.sqanti3_chunk_manifest.tsv"), emit: sqanti3_results

    when:
        task.ext.when == null || task.ext.when

    script:
        def chunk_suffix = chunk_id == "full" ? "" : "_${chunk_id}"
        def output_prefix = "${params.runID}_${genome_id}_${stage}${chunk_suffix}"
        def include_orf = params.sqanti3_include_orf ? "--include_ORF" : ""
        def tusco = params.sqanti3_tusco ? "--tusco ${params.sqanti3_tusco}" : ""
        def polya_available = polya_source && !polya_source.toString().toLowerCase().matches(/^(|na|none|null|\.)$/)
        def cage_available = cage_source && !cage_source.toString().toLowerCase().matches(/^(|na|none|null|\.)$/)
        def polya = polya_available ? "--polyA_peak ${polya_peak}" : ""
        def cage = cage_available ? "--CAGE_peak ${cage_peak}" : ""
        def sqanti3_by_chromosome = params.sqanti3_by_chromosome.toString().toLowerCase() in ['true', '1', 'yes']
        def report = sqanti3_by_chromosome ? "--report skip" : "--report ${params.sqanti3_report}"
        def output_dir = "${stage}_${chunk_id}"
        """
        sqanti_out="\$PWD/sqanti3_output/chunks/${output_dir}"
        mkdir -p "\$sqanti_out"
        manifest="\$(mktemp ${output_prefix}.XXXXXX.sqanti3_chunk_manifest.tsv)"

        conda run -n sqanti3 sqanti3_qc.py \\
            --isoforms ${isoforms} \\
            --refGTF ${reference_gtf} \\
            --refFasta ${reference_fasta} \\
            ${include_orf} \\
            ${tusco} \\
            ${polya} \\
            ${cage} \\
            -o ${output_prefix} \\
            -d "\$sqanti_out" \\
            --cpus ${task.cpus} \\
            ${report}

        {
            printf 'stage\\tgenome_id\\tchunk_id\\toutput_dir\\n'
            printf '${stage}\\t${genome_id}\\t${chunk_id}\\t%s\\n' "\$sqanti_out"
        } >| "\$manifest"
        """
}
