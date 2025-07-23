//
// Subworkflow for the extraction of transcript models from bam files for the LoRA pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAM_TO_BED      } from "../../../modules/local/bam_to_bed/main.nf"
include { BED_STRAND_FIX  } from '../../../modules/local/bed_strand_fix/main.nf'
include { BED_TO_GENEPRED } from '../../../modules/local/bed_to_genepred/main.nf'
include { GENEPRED_TO_GTF } from '../../../modules/local/genepred_to_gtf/main.nf'
include { FIX_GTF         } from '../../../modules/local/fix_gtf/main.nf'
include { FILTER_BAM_FOR_CANONICAL } from '../../../modules/local/filter_bam_for_canonical/main.nf'
include { SAMTOOLS_MERGE_CANONICAL } from '../../../modules/local/samtools_merge_canonical/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BAM_TO_GTF {
    take:
    bam_files
    chromosome_bam

    main:
    // Converting bam files to bed format
    BAM_TO_BED(bam_files)

    // Attempting to fix strand for multi-exon transcripts
    BED_STRAND_FIX(BAM_TO_BED.out.converted_bed)

    // Filter bam file to extract only canonical splicing reads
    FILTER_BAM_FOR_CANONICAL(chromosome_bam.map { bam, _bai, _reference, genome_id, chromosome_id, chromosome_count, _bin_count -> tuple(chromosome_id, genome_id, chromosome_count, bam) }.join(BED_STRAND_FIX.out.fixed_bed.map { bed, genome_id, chromosome_id, chromosome_count, bin_count -> tuple(groupKey(chromosome_id, bin_count.toInteger()), bed, genome_id, chromosome_count) }.groupTuple().map { chromosome_id, bed, genome_id, chromosome_count -> [chromosome_id, genome_id[0], chromosome_count[0], bed] }, by: [0,1,2]))

    // Merge the filtered bam files
    SAMTOOLS_MERGE_CANONICAL(FILTER_BAM_FOR_CANONICAL.out.chromosome_canonical_bam.map { _chromosome_id, genome_id, chromosome_count, canonical_bam, canonical_bam_bai -> tuple(groupKey(genome_id, chromosome_count.toInteger()), canonical_bam, canonical_bam_bai) }.groupTuple())

    // Converting bed format to genepred format
    BED_TO_GENEPRED(BED_STRAND_FIX.out.fixed_bed)

    // Converting genePred format to GTF format
    GENEPRED_TO_GTF(BED_TO_GENEPRED.out.converted_genepred)

    // Prepare GTF for first round of TMERGE
    FIX_GTF(GENEPRED_TO_GTF.out.converted_gtf)

    emit:
    prepared_gtf = FIX_GTF.out.prepared_gtf
    chromosome_bam = FILTER_BAM_FOR_CANONICAL.out.chromosome_canonical_bam.map { chromosome_id, genome_id, chromosome_count, canonical_bam, canonical_bam_bai -> [canonical_bam, canonical_bam_bai, genome_id, chromosome_id, chromosome_count] }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/