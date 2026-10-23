process PREPARE_SQANTI3_REFERENCE {
    label "process_low"
    label "error_retry"

    input:
        tuple path(reference_fasta), val(genome_id), path(reference_gtf), val(polya_peak), val(cage_peak)

    output:
        tuple val(genome_id),
            path("${genome_id}.sqanti3.reference.fa"),
            path("${genome_id}.sqanti3.reference.gtf"),
            path("${genome_id}.sqanti3.polya_peak.bed"),
            val(polya_peak),
            path("${genome_id}.sqanti3.cage_peak.bed"),
            val(cage_peak),
            emit: prepared_reference

    when:
        task.ext.when == null || task.ext.when

    script:
        def polya_available = polya_peak && !polya_peak.toString().toLowerCase().matches(/^(|na|none|null|\.)$/)
        def cage_available = cage_peak && !cage_peak.toString().toLowerCase().matches(/^(|na|none|null|\.)$/)
        """
        if [[ "${reference_fasta}" == *.gz ]]; then
            gzip -dc "${reference_fasta}" > ${genome_id}.sqanti3.reference.fa
        else
            ln -s "${reference_fasta}" ${genome_id}.sqanti3.reference.fa
        fi

        if [[ "${reference_gtf}" == *.gz ]]; then
            gzip -dc "${reference_gtf}" > ${genome_id}.sqanti3.reference.gtf
        else
            ln -s "${reference_gtf}" ${genome_id}.sqanti3.reference.gtf
        fi

        if [[ "${polya_available}" == "true" ]]; then
            if [[ "${polya_peak}" == *.gz ]]; then
                gzip -dc "${polya_peak}" > ${genome_id}.sqanti3.polya_peak.bed
            else
                ln -s "${polya_peak}" ${genome_id}.sqanti3.polya_peak.bed
            fi
        else
            touch ${genome_id}.sqanti3.polya_peak.bed
        fi

        if [[ "${cage_available}" == "true" ]]; then
            if [[ "${cage_peak}" == *.gz ]]; then
                gzip -dc "${cage_peak}" > ${genome_id}.sqanti3.cage_peak.bed
            else
                ln -s "${cage_peak}" ${genome_id}.sqanti3.cage_peak.bed
            fi
        else
            touch ${genome_id}.sqanti3.cage_peak.bed
        fi
        """
}
