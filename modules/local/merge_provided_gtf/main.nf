process MERGE_PROVIDED_GTF {
    label 'process_very_low'
    label 'error_retry'
    input:
        tuple path(input_gtfs), val(genome_id), val(individual_bin_count), val(total_bin_count)
    output:
        tuple path("${params.runID}_fixed.gtf"), val(genome_id), val(individual_bin_count), val(total_bin_count), emit: merged_fixed_gtf
    script:
        """
        cat ${input_gtfs} > ${params.runID}.gtf
        awk -F'\t' '\$5-\$4 > 0 && \$3 == "exon" {{print}}' ${params.runID}.gtf > ${params.runID}_fixed.gtf
        """
}
