process SPLIT_GTF_BY_CHROMOSOME {
    label 'process_very_low'
    label 'error_retry'
    input:
        tuple val(chromosome), val(genome_id), path(tmerge_gtf), path(contains_file), val(chromosome_count)
    output:
        tuple val(chromosome), val(genome_id), path("${tmerge_gtf.simpleName}_${chromosome}.gtf"), path("${contains_file.simpleName}_${chromosome}.tsv"), val(chromosome_count), emit: chrom_gtf
    script:
        """
        chromosome_id=${chromosome}
        grep -E "^\\b\$chromosome_id\\b" ${tmerge_gtf} > ${tmerge_gtf.simpleName}_\$chromosome_id.gtf
        awk -v chromosome_id="\$chromosome_id" '\$1 == chromosome_id {print \$0}' ${contains_file} > ${contains_file.simpleName}_\$chromosome_id.tsv
        """
}
