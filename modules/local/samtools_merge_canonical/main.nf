process SAMTOOLS_MERGE_CANONICAL {
    label "process_medium"
    label "multi_cpu_high"
    label "error_retry"
    container "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    publishDir "${params.outDir}/${genome_id}/reads/canonical/merged/", mode: 'copy', pattern: "*.canonical.bam*"
    input:
        tuple val(genome_id), path(bams), path(bai)
    output:
        tuple path("${params.runID}_${genome_id}.canonical.bam"), path("${params.runID}_${genome_id}.canonical.bam.csi"), val(genome_id), emit: merged_bam
    script:
        """
        samtools merge -@ 20 -o ${params.runID}_${genome_id}.canonical.bam ${bams}
        samtools index -c -@ 20 ${params.runID}_${genome_id}.canonical.bam
        """
}