process EXTRACT_CHROMOSOME_IDS_FROM_GTF {
    label 'process_very_low'
    label 'error_retry'
    input:
    tuple val(genome_id), path(tmerged_gtf), path(contains_file)
    output:
    path('tmerge_sample_sheet.csv'), emit: tmerge_sample_sheet
    when:
        task.ext.when == null || task.ext.when

    script:
    """
    source ~/venv/bin/activate
    extract_chromosomes_from_gtf.py -g ${genome_id} -f ${tmerged_gtf} -c ${contains_file} -o tmerge_sample_sheet.csv
    """
}