process FIX_GTF {
    label "process_low"
    label "error_retry"

    input:
    tuple path(gtf), val(genome_id), val(chromosome_id), val(chromosome_count), val(bin_count)

    output:
    tuple path("*.fixed.gtf"), val(genome_id), val(chromosome_id), val(chromosome_count), val(bin_count), emit: prepared_gtf

    script:
    """
    awk -F'\t' '\$5-\$4 > 0 && \$3 == "exon" {print}' ${gtf} > temp.gtf
    sort -V -k1,1 -k4,4n temp.gtf > "${gtf.simpleName}_bam.fixed.gtf"
    """
}