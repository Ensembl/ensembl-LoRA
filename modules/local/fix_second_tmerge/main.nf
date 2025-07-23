process FIX_SECOND_TMERGE {
    label 'process_very_low'
    label 'error_retry'
    input:
        tuple val(genome_id), path(input_gtf), path(initial_merged_contains), val(chromosome_count), val(chromosome_id)
    output:
        tuple val(genome_id), path("*.retmerge.gtf"), path(initial_merged_contains), val(chromosome_count), val(chromosome_id), emit: second_tmerge_fixed
    script:
        """
        awk 'BEGIN{OFS=FS="\\t"; IGNORECASE=1} {gsub(/TM_/, "chr"\$1"_", \$NF); gsub(/chrchr/, "chr", \$NF); print}' ${input_gtf} > ${input_gtf.simpleName}.retmerge.gtf
        """
}
