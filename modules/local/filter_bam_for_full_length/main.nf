process FILTER_BAM_FOR_FULL_LENGTH {
    label 'process_low'
    label 'multi_cpu_high'
    label 'error_retry'
    publishDir "${params.outDir}/${genome_id}/full_length/reads/chromosomes/", mode: 'copy', pattern: "*.full_length.bam*"
    input:
        tuple val(chromosome_id), val(genome_id), val(chromosome_count), path(bam), path(full_length_contains)
    output:
        tuple val(genome_id), val(chromosome_count), path("${bam.simpleName}.full_length.bam"), path(full_length_contains), emit: full_length_filtered_bam
        tuple val(genome_id), val(chromosome_count), path("${bam.simpleName}.full_length.bam.csi"), emit: full_length_filtered_index
    script:
    """
    cut -f3 ${full_length_contains} | tr ',' '\\n' > reads_list.txt
    samtools view -H ${bam} > header.sam
    grep -v "^@PG" header.sam > header.sam.tmp
    mv header.sam.tmp header.sam
    if [ -s reads_list.txt ]; then
        samtools view -@ 20 ${bam} > temp.sam
        LC_ALL=C grep -Fwf reads_list.txt temp.sam > filtered_reads.sam
        cat header.sam filtered_reads.sam | samtools view -b -@ 20 > ${bam.simpleName}.full_length.bam
        samtools index -c -@ 20 ${bam.simpleName}.full_length.bam
    else
        samtools view -bS header.sam > ${bam.simpleName}.full_length.bam
        samtools index -c -@ 20 ${bam.simpleName}.full_length.bam
    fi
    rm -f reads_list.txt header.sam temp.sam filtered_reads.sam
    """
}
// samtools sort -@ 20 -O bam -o ${bam.simpleName}.full_length.bam ${bam.simpleName}.temp.bam
// container "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"