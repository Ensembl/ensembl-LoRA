process FILTER_BAM_FOR_FULL_LENGTH {
    label 'process_low'
    label 'multi_cpu_high'
    label 'error_retry'
    input:
        tuple val(chromosome_id), val(genome_id), val(chromosome_count), path(bam), path(full_length_contains)
    output:
        tuple val(genome_id), val(chromosome_count), path("${bam.simpleName}.full_length.bam"), path(full_length_contains), emit: full_length_filtered_bam
        tuple val(genome_id), val(chromosome_count), path("${bam.simpleName}.full_length.bam.csi"), emit: full_length_filtered_index
    when:
        task.ext.when == null || task.ext.when

    script:
    """
    cut -f3 ${full_length_contains} | tr ',' '\\n' > reads_list.txt
    if [ -s reads_list.txt ]; then
        samtools view -@ ${task.cpus} -N reads_list.txt -b -o ${bam.simpleName}.full_length.bam ${bam}
        samtools index -c -@ ${task.cpus} ${bam.simpleName}.full_length.bam
    else
        samtools view -@ ${task.cpus} -H ${bam} | grep -v "^@PG" | samtools view -@ ${task.cpus} -b -o ${bam.simpleName}.full_length.bam -
        samtools index -c -@ ${task.cpus} ${bam.simpleName}.full_length.bam
    fi
    rm -f reads_list.txt
    """
}
// samtools sort -@ 20 -O bam -o ${bam.simpleName}.full_length.bam ${bam.simpleName}.temp.bam
// container "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
