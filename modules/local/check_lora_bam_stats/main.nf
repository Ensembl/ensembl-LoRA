process CHECK_LORA_BAM_STATS {
    label "process_medium"
    label "multi_cpu_low"
    label "error_retry"
    container "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"

    input:
        tuple val(stage), val(genome_id), val(label), path(bam), path(bai)

    output:
        tuple val(genome_id), path("${params.runID}_${genome_id}_${stage}_${label}_bam_stats.tsv"), emit: bam_stats

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        samtools flagstat -@ ${task.cpus} -O tsv ${bam} > flagstat.tsv
        total=\$(awk -F '\\t' '\$3 == "total" { print \$1; exit }' flagstat.tsv)
        mapped=\$(awk -F '\\t' '\$3 == "mapped" { print \$1; exit }' flagstat.tsv)
        secondary=\$(awk -F '\\t' '\$3 == "secondary" { print \$1; exit }' flagstat.tsv)
        supplementary=\$(awk -F '\\t' '\$3 == "supplementary" { print \$1; exit }' flagstat.tsv)
        total=\${total:-0}
        mapped=\${mapped:-0}
        secondary=\${secondary:-0}
        supplementary=\${supplementary:-0}
        unmapped=\$(( total - mapped ))
        low_mapq=\$(samtools view -@ ${task.cpus} ${bam} | awk '\$5 < 30 { count += 1 } END { print count + 0 }')
        {
            printf "stage\\tgenome_id\\tlabel\\tbam\\tindex\\ttotal\\tmapped\\tunmapped\\tsecondary\\tsupplementary\\tlow_mapq\\tbam_size_bytes\\tindex_size_bytes\\n"
            printf "${stage}\\t${genome_id}\\t${label}\\t${bam}\\t${bai}\\t\$total\\t\$mapped\\t\$unmapped\\t\$secondary\\t\$supplementary\\t\$low_mapq\\t%s\\t%s\\n" "\$(stat -c %s ${bam})" "\$(stat -c %s ${bai})"
        } > ${params.runID}_${genome_id}_${stage}_${label}_bam_stats.tsv
        """
}
