process PB_MINIMAP2 {
    label "process_very_high"
    label "multi_cpu_medium"
    label "error_retry"
    container "docker://community.wave.seqera.io/library/minimap2_samtools:3a608f92b72cacce"
    input:
        tuple path(reads), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species)
    output:
        tuple path("*.sorted.bam"), path("*.sorted.bam.csi"), path(reference), val(genome_id), val(raw_file_count), emit: mapped_reads_pb
    when:
        task.ext.when == null || task.ext.when

    script:
        """
        if [ ${params.genomeSize} == "normal" ]; then
            minimap2 --MD -L -ax splice:hq --secondary=no -t ${task.cpus} ${reference} ${reads} | samtools view -@ ${task.cpus} -F 256 -F 4 -F 2048 -q 30 -b - | samtools sort -@ ${task.cpus} -O bam -o ${reads.simpleName}_${genome_id}.sorted.bam -
        else
            minimap2 --MD --split-prefix=tmp -L -ax splice:hq --secondary=no -t ${task.cpus} ${reference} ${reads} | samtools view -@ ${task.cpus} -F 256 -F 4 -F 2048 -q 30 -b - | samtools sort -@ ${task.cpus} -O bam -o ${reads.simpleName}_${genome_id}.sorted.bam -
        fi
        samtools index -c -@ ${task.cpus} ${reads.simpleName}_${genome_id}.sorted.bam
        """
}
