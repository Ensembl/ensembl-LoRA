process FIX_GTF {
    label "process_low"
    label "error_retry"

    input:
    tuple path(gtf), val(genome_id), val(chromosome_id), val(chromosome_count), val(bin_count)

    output:
    tuple path("*.fixed.gtf"), val(genome_id), val(chromosome_id), val(chromosome_count), val(bin_count), emit: prepared_gtf

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        awk -F'\t' '\$5-\$4 > 0 && \$3 == "exon" {print}' ${gtf} > temp.gtf
        sort -V -k1,1 -k4,4n temp.gtf > "${gtf.simpleName}_bam.fixed.gtf"
        if [ \$? -eq 0 ]; then
            rm -f temp.gtf
            rm -f ${gtf} \$(readlink -f ${gtf})
        fi
        """
}