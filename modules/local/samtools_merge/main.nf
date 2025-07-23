process SAMTOOLS_MERGE {
    label "process_medium"
    label "multi_cpu_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    publishDir "${params.outDir}/${genome_id}/reads/all_reads/merged/", mode: 'copy', pattern: "*.bam*"
    input:
        tuple path(bams), path(reference_genome), val(genome_id)
    output:
        tuple path("${params.runID}_${genome_id}.bam"), path("${params.runID}_${genome_id}.bam.csi"), path(reference_genome), val(genome_id), emit: merged_bam
    script:
        """
        samtools merge -@ 20 -o ${params.runID}_${genome_id}.bam ${bams}
        samtools index -c -@ 20 ${params.runID}_${genome_id}.bam
        """
}