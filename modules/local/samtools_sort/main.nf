process SAMTOOLS_SORT {
    label "process_high"
    label "multi_cpu_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    publishDir "${params.outDir}/${genome_id}/reads/all_reads/split/", mode: 'copy', pattern: "*.bam*"
    input:
        tuple path(bam), path(reference_genome), val(genome_id), val(raw_file_count)
    output: 
        tuple val(genome_id), path("*.bam"), path(reference_genome), val(raw_file_count), emit: sorted_bam
        tuple val(genome_id), path("*.bam.csi"), path(reference_genome), val(raw_file_count), emit: sorted_bam_bai
    script:
        """
        samtools sort -@ 20 -O bam -o ${bam.simpleName}_${genome_id}.sorted.bam ${bam}
        samtools index -c -@ 20 ${bam.simpleName}_${genome_id}.sorted.bam
        """
}