process FIX_FIRST_TMERGE {
    label 'process_very_low'
    label 'error_retry'
    input:
        tuple val(chromosome_id), path(input_gtf), val(genome_id), val(chromosome_count), val(bin_count)
    output:
        tuple val(chromosome_id), path("*.tmerge.gtf"), val(genome_id), val(chromosome_count), val(bin_count), emit: first_tmerge_fixed
    script:
        """
        awk 'BEGIN{OFS=FS="\\t"; IGNORECASE=1} {gsub(/TM_/, "chr"\$1"_", \$NF); gsub(/chrchr/, "chr", \$NF); print}' ${input_gtf} > ${input_gtf.simpleName}.tmerge.gtf
        """
}
