process COMBINE_READ_ID_MAPS {
    label 'process_medium'
    label 'error_retry'
    input:
        tuple val(genome_id), path(read_id_maps)
    output:
        tuple val(genome_id), path("${params.runID}_${genome_id}_read_id_map.tsv.gz"), emit: combined_read_id_map
    when:
        task.ext.when == null || task.ext.when

    script:
    """
    first=1
    for read_map in *.read_id_map.tsv *.read_id_map.tsv.gz; do
        [ -e "\$read_map" ] || continue
        if [ "\$first" -eq 1 ]; then
            gzip -cdf "\$read_map"
            first=0
        else
            gzip -cdf "\$read_map" | awk 'NR > 1'
        fi
    done | gzip -c > ${params.runID}_${genome_id}_read_id_map.tsv.gz
    """
}
