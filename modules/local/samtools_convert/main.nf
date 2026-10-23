process SAMTOOLS_CONVERT {
    label "process_low"
    label "multi_cpu_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    input:
        tuple path(aligned_reads), path(reference_genome), val(genome_id), val(raw_file_count)
    output:
        tuple path("*.bam"), path(reference_genome), val(genome_id), val(raw_file_count), emit: converted_bam
    when:
        task.ext.when == null || task.ext.when

    script:
        """
        samtools view -@ 20 -S -b -o temp ${aligned_reads}
        samtools view -@ 20 -F 256 -F 4 -F 2048 -q 30 -bo ${aligned_reads.simpleName}.bam temp
        if [ \$? -eq 0 ]; then
            rm -f temp
            rm -f ${aligned_reads} \$(readlink -f ${aligned_reads})
        fi
        """
}
