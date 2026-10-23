process BED_STRAND_FIX {
    label "process_very_low"
    label "error_retry"
    container "docker://quay.io/biocontainers/perl-bio-db-hts:3.01--pl5321h6141fd1_8"
    input:
        tuple path(bed), path(reference_genome), val(genome_id), val(chromosome_id), val(chromosome_count), val(bin_count)
    output:
        tuple path("*fixed*.bed"), val(genome_id), val(chromosome_id), val(chromosome_count), val(bin_count), emit: fixed_bed
        path("*log"), emit: log
    when:
        task.ext.when == null || task.ext.when

    script:
        """
        source ~/venv/bin/activate
        touch ${bed.simpleName}.fixed.bed
        splice_site_checker.py ${bed} ${reference_genome} ${bed.simpleName}
        if [ ! -s ${bed.simpleName}.fixed.bed ]; then
            touch ${bed.simpleName}.fixed.bed
        fi
        rm -f ${reference_genome}.*
        """
}
// source ~/venv/bin/activate
// gunzip -c ${reference_genome} > reference.fa
// touch ${bed.simpleName}.fixed.bed
// splice_site_checker.py ${bed} reference.fa ${bed.simpleName}
// if [ ! -s ${bed.simpleName}.fixed.bed ]; then
//     touch ${bed.simpleName}.fixed.bed
// fi