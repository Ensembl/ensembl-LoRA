//
// Subworkflow for the preparation of mapped reads for the LoRA pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_SORT        } from '../../../modules/local/samtools_sort/main.nf'
include { SAMTOOLS_MERGE       } from '../../../modules/local/samtools_merge/main.nf'
include { PICARD_SPLIT         } from '../../../modules/local/picard_split/main.nf'
include { SPLIT_PICARD_CHANNEL } from '../../../modules/local/split_picard_channel/main.nf'
include { CHROMOSOME_SAMPLE_SHEET } from '../../../modules/local/chromosome_sample_sheet/main.nf'
include { CHROMOSOME_SPLIT } from '../../../modules/local/chromosome_split/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SAMTOOLS {
    take:
    mapped_bam_files
    bam_files

    main:
    // Minimap2-derived BAMs are already sorted and indexed. User-provided BAMs still need normalisation.
    SAMTOOLS_SORT(bam_files)

    sorted_mapped_bam = mapped_bam_files
        .map { bam, bai, reference, genome_id, raw_file_count ->
            tuple(genome_id, bam, bai, reference, raw_file_count)
        }

    sorted_bam_with_index = sorted_mapped_bam.mix(SAMTOOLS_SORT.out.sorted_bam_with_index)

    sorted_bam_for_merge = sorted_bam_with_index
        .map { genome_id, bam, _bai, reference, raw_file_count ->
            tuple(genome_id, bam, reference, raw_file_count)
        }

    // Merge all bam files together
    SAMTOOLS_MERGE(sorted_bam_for_merge.map { genome_id, bam, reference, raw_file_count -> tuple(groupKey(genome_id, raw_file_count.toInteger()), bam, reference)}.groupTuple().map { genome_id, bams, fasta -> [bams, fasta[0], genome_id] })

    // Creating chromosome sample sheet
    CHROMOSOME_SAMPLE_SHEET(SAMTOOLS_MERGE.out.merged_bam)

    // Split each bam into chroosomes
    CHROMOSOME_SPLIT(CHROMOSOME_SAMPLE_SHEET.out.chromosome_sample_sheet.splitCsv())

    // Split each chromosome bam into bins
    PICARD_SPLIT(CHROMOSOME_SPLIT.out.chromosome_bam)

    // Split PICARD_SPLIT output into separate channels
    SPLIT_PICARD_CHANNEL(PICARD_SPLIT.out.picard_split)

    emit:
    prepared_bam = SPLIT_PICARD_CHANNEL.out.picard_sample_sheet.splitCsv()
    chromosome_bam = CHROMOSOME_SPLIT.out.chromosome_bam
    sorted_bam = sorted_bam_for_merge
    sorted_bam_bai = sorted_bam_with_index.map { genome_id, _bam, bai, reference, raw_file_count -> tuple(genome_id, bai, reference, raw_file_count) }
    sorted_bam_with_index = sorted_bam_with_index
    merged_bam = SAMTOOLS_MERGE.out.merged_bam
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
