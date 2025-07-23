process CHROMOSOME_SAMPLE_SHEET {
    label "process_very_low"
    label "multi_cpu_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    input:
        tuple path(bam), path(bai), path(reference_genome), val(genome_id)
    output:
        path("${bam.simpleName}_sample_sheet.csv"), emit: chromosome_sample_sheet
    script:
        """
        source ~/venv/bin/activate
        chromosome_sample_sheet.py ${genome_id} ${bam} ${bai} ${reference_genome} 100000 ${bam.simpleName}
        """
}