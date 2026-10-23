//
// Subworkflow for recreating split BAM channels from a published merged BAM
//

include { CHROMOSOME_SAMPLE_SHEET } from '../../../modules/local/chromosome_sample_sheet/main.nf'
include { CHROMOSOME_SPLIT        } from '../../../modules/local/chromosome_split/main.nf'
include { PICARD_SPLIT            } from '../../../modules/local/picard_split/main.nf'
include { SPLIT_PICARD_CHANNEL    } from '../../../modules/local/split_picard_channel/main.nf'

workflow SAMTOOLS_FROM_MERGED_BAM {
    take:
    merged_bam

    main:
    CHROMOSOME_SAMPLE_SHEET(merged_bam)
    CHROMOSOME_SPLIT(CHROMOSOME_SAMPLE_SHEET.out.chromosome_sample_sheet.splitCsv())
    PICARD_SPLIT(CHROMOSOME_SPLIT.out.chromosome_bam)
    SPLIT_PICARD_CHANNEL(PICARD_SPLIT.out.picard_split)

    emit:
    prepared_bam = SPLIT_PICARD_CHANNEL.out.picard_sample_sheet.splitCsv()
    chromosome_bam = CHROMOSOME_SPLIT.out.chromosome_bam
}
