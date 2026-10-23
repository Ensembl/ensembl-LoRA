process FIRST_TMERGE {
    label 'process_very_low'
    label 'error_retry'
    input:
        tuple path(input_gtf), val(genome_id), val(chromosome_id), val(chromosome_count), val(bin_count)
    output:
        tuple val(chromosome_id), path("*.unfixed.gtf"), val(genome_id), val(chromosome_count), val(bin_count), emit: first_tmerge
    when:
        task.ext.when == null || task.ext.when

    script:
        """
        tmerge.pl --exonOverhangTolerance 10 ${input_gtf} > ${input_gtf.simpleName}.unfixed.gtf
        """
}