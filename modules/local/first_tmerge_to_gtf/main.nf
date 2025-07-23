process FIRST_TMERGE_TO_GTF {
    label "process_very_low"
    label "error_retry"
    input:
        tuple val(chromosome_id), path(tmerge_gtf_ch), val(genome_id), val(chromosome_count), val(bin_count)
    output:
        tuple val(chromosome_id), path("*.tama.gtf"), path("*.contains.tsv"), val(genome_id), val(chromosome_count), val(bin_count), emit: first_tmerge_gtf
    script:
        """
        if [ -s ${tmerge_gtf_ch} ]; then
            source ~/venv/bin/activate
            tmerge_to_gtf.py ${tmerge_gtf_ch} ${tmerge_gtf_ch.simpleName} ${tmerge_gtf_ch.simpleName}
        else
            touch ${tmerge_gtf_ch.simpleName}.tama.gtf
            touch ${tmerge_gtf_ch.simpleName}.contains.tsv
        fi
        """
}
