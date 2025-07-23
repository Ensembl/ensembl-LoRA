process PB_MINIMAP2 {
    label "process_very_high"
    label "multi_cpu_medium"
    label "error_retry"
    container "docker://quay.io/biocontainers/minimap2:2.26--he4a0461_2"
    input:
        tuple path(reads), val(seq_source), path(reference), val(genome_id), val(raw_file_count), val(species)
    output:
        tuple path("*.sam"), path(reference), val(genome_id), val(raw_file_count), emit: mapped_reads_pb
    script:
        """
        if [ ${species} == "human" ]; then
            minimap2 --MD -L -ax splice:hq --secondary=no -t 10 ${reference} ${reads} > ${reads.simpleName}.sam
        else
            minimap2 --split-prefix=tmp -L -ax splice:hq --secondary=no -t 10 ${reference} ${reads} > ${reads.simpleName}.sam
        fi
        """
}   