process PICARD_SPLIT {
    label "process_medium"
    label "error_retry"
    input:
        tuple path(bam), path(bam_bai), path(reference), val(genome_id), val(chromosome_id), val(chromosome_count), val(bin_count)
    output:
        tuple path("*.picard.bam"), path(reference), val(genome_id), val(chromosome_id), val(chromosome_count), val(bin_count), emit: picard_split
    script:
        """
        if [ ${bin_count} -gt 1 ]; then
            picard SplitSamByNumberOfReads --INPUT ${bam} --OUTPUT ./ --SPLIT_TO_N_FILES ${bin_count}
            source ~/venv/bin/activate
            rename_shards.py -s shard* -b ${bam}
        else
            cp ${bam} ${bam.simpleName}_0001.picard.bam
        fi
        """
}
