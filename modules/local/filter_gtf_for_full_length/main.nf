process FILTER_GTF_FOR_FULL_LENGTH {
    label 'process_very_low'
    label 'error_retry'
    input:
        tuple val(chromosome_id), val(genome_id), val(chromosome_count), path(gtf), path(full_length_contains)
    output:
        tuple val(genome_id), val(chromosome_count), val(chromosome_id), path("${gtf.simpleName}.filtered.gtf"), path("${params.runID}_${chromosome_id}_FL_contains.tsv"), emit: full_length_filtered_gtf
    when:
        task.ext.when == null || task.ext.when

    script:
    """
    awk -F',' 'NF > x' x=\$((${params.minFL} - 1)) ${full_length_contains} > ${params.runID}_${chromosome_id}_FL_contains.tsv
    if [ -s ${params.runID}_${chromosome_id}_FL_contains.tsv ]; then
        cut -f2 ${params.runID}_${chromosome_id}_FL_contains.tsv > patterns.txt
    grep -Fwf patterns.txt ${gtf} > ${gtf.simpleName}.filtered.gtf
    else
        touch ${gtf.simpleName}.filtered.gtf
    fi
    """
}
