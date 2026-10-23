process PREPARE_SQANTI3_CHROMOSOME_REFERENCE {
    label "process_low"
    label "error_retry"

    input:
        tuple val(genome_id),
            val(stage),
            val(chromosome_id),
            path(isoforms),
            path(reference_fasta),
            path(reference_gtf),
            path(polya_peak),
            val(polya_source),
            path(cage_peak),
            val(cage_source)

    output:
        tuple val(stage),
            val(genome_id),
            val(chromosome_id),
            path(isoforms),
            path("${genome_id}_${chromosome_id}.sqanti3.reference.fa"),
            path("${genome_id}_${chromosome_id}.sqanti3.reference.gtf"),
            path("${genome_id}_${chromosome_id}.sqanti3.polya_peak.bed"),
            val(polya_source),
            path("${genome_id}_${chromosome_id}.sqanti3.cage_peak.bed"),
            val(cage_source),
            path("${genome_id}_${chromosome_id}.sqanti3.chromosome_map.tsv"),
            emit: chromosome_reference

    when:
        task.ext.when == null || task.ext.when

    script:
        def polya_available = polya_source && !polya_source.toString().toLowerCase().matches(/^(|na|none|null|\.)$/)
        def cage_available = cage_source && !cage_source.toString().toLowerCase().matches(/^(|na|none|null|\.)$/)
        def polya_filter = polya_available ?
            "awk -v chrom=\"\$sqanti_chrom\" 'BEGIN { FS=OFS=\"\\t\" } \$1 == chrom { print }' ${polya_peak} > ${genome_id}_${chromosome_id}.sqanti3.polya_peak.bed" :
            "touch ${genome_id}_${chromosome_id}.sqanti3.polya_peak.bed"
        def cage_filter = cage_available ?
            "awk -v chrom=\"\$sqanti_chrom\" 'BEGIN { FS=OFS=\"\\t\" } \$1 == chrom { print }' ${cage_peak} > ${genome_id}_${chromosome_id}.sqanti3.cage_peak.bed" :
            "touch ${genome_id}_${chromosome_id}.sqanti3.cage_peak.bed"
        """
        ln -s ${reference_fasta} ${genome_id}_${chromosome_id}.sqanti3.reference.fa

        sqanti_chrom=\$(awk 'BEGIN { FS="\\t" } \$0 !~ /^#/ && NF >= 9 { print \$1; exit }' ${isoforms})
        if [ -z "\$sqanti_chrom" ]; then
            echo "No feature rows found in SQANTI3 isoform GTF for LoRA chromosome_id '${chromosome_id}': ${isoforms}" >&2
            exit 1
        fi

        chromosome_count=\$(awk 'BEGIN { FS="\\t" } \$0 !~ /^#/ && NF >= 9 { seen[\$1] = 1 } END { for (chrom in seen) count += 1; print count + 0 }' ${isoforms})
        if [ "\$chromosome_count" -ne 1 ]; then
            echo "SQANTI3 chromosome-mode expected one real chromosome in ${isoforms} for LoRA chromosome_id '${chromosome_id}', but found \$chromosome_count." >&2
            awk 'BEGIN { FS="\\t" } \$0 !~ /^#/ && NF >= 9 { seen[\$1] = 1 } END { for (chrom in seen) print chrom }' ${isoforms} >&2
            exit 1
        fi

        {
            echo -e "chromosome_id\\tsqanti_chrom"
            echo -e "${chromosome_id}\\t\$sqanti_chrom"
        } > ${genome_id}_${chromosome_id}.sqanti3.chromosome_map.tsv

        awk -v chrom="\$sqanti_chrom" 'BEGIN { FS=OFS="\\t" } /^#/ { print; next } \$1 == chrom { print }' ${reference_gtf} > ${genome_id}_${chromosome_id}.sqanti3.reference.gtf
        reference_feature_count=\$(awk 'BEGIN { FS="\\t" } \$0 !~ /^#/ && NF >= 9 { count += 1 } END { print count + 0 }' ${genome_id}_${chromosome_id}.sqanti3.reference.gtf)
        if [ "\$reference_feature_count" -eq 0 ]; then
            echo "Reference GTF ${reference_gtf} does not contain records for model chromosome '\$sqanti_chrom' from LoRA chromosome_id '${chromosome_id}'." >&2
            echo "Check that LoRA model GTF chromosome names match the SQANTI3 reference GTF/BED chromosome names." >&2
            exit 1
        fi

        ${polya_filter}
        ${cage_filter}
        """
}
