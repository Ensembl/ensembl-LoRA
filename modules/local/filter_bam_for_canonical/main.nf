process FILTER_BAM_FOR_CANONICAL {
    label "process_very_low"
    label "multi_cpu_low"
    label "error_retry"
    input:
        tuple val(chromosome_id), val(genome_id), val(chromosome_count), path(bam), path(bed)
    output:
        tuple val(chromosome_id), val(genome_id), val(chromosome_count), path("${params.runID}_${chromosome_id}.canonical.bam"), path("${params.runID}_${chromosome_id}.canonical.bam.csi"), emit: chromosome_canonical_bam
    when:
        task.ext.when == null || task.ext.when

    script:
        """
        source ~/venv/bin/activate
        filter_and_fix_bams.py --bam ${bam} --prefix ${params.runID}_${chromosome_id} --threads ${task.cpus}
        if [ \$? -eq 0 ]; then
            rm -f ${bam} \$(readlink -f ${bam})
        fi
        """
}
