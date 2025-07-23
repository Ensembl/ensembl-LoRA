//
// Subworkflow for the merging of transcripts for the LoRA pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INTRON_CHAIN_EXTRACTOR_READS } from '../../../modules/local/intron_chain_extractor_reads/main.nf'
include { INTRON_CHAIN_EXTRACTOR_MODELS } from '../../../modules/local/intron_chain_extractor_models/main.nf'
include { COMBINE_INTRON_CHAINS_READS } from '../../../modules/local/combine_intron_chains_reads/main.nf'
include { COMBINE_INTRON_CHAINS_MODELS } from '../../../modules/local/combine_intron_chains_models/main.nf'
include { FULL_LENGTH_CONTAINS } from '../../../modules/local/full_length_contains/main.nf'
include { BAM_TO_BED_FULL_LENGTH_FILTER } from '../../../modules/local/bam_to_bed_full_length_filter/main.nf'
include { FILTER_GTF_FOR_FULL_LENGTH } from '../../../modules/local/filter_gtf_for_full_length/main.nf'
include { FILTER_BAM_FOR_FULL_LENGTH } from '../../../modules/local/filter_bam_for_full_length/main.nf'
include { COMBINE_FILTERED_MODELS } from '../../../modules/local/combine_filtered_models/main.nf'
include { COMBINE_FULL_LENGTH_READS } from '../../../modules/local/combine_full_length_reads/main.nf'
include { COMBINE_FULL_LENGTH_STATS } from '../../../modules/local/combine_full_length_stats/main.nf'
include { COMBINE_UNFILTERED_FULL_LENGTH_CONTAINS } from '../../../modules/local/combine_unfiltered_full_length_contains/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FULL_LENGTH_FILTER {
    take:
    merged_reads
    merged_transcripts

    main:
    // Convert bam to bed
    BAM_TO_BED_FULL_LENGTH_FILTER(merged_reads)

    // Extract intron chains and single exons from merged reads
    INTRON_CHAIN_EXTRACTOR_READS(BAM_TO_BED_FULL_LENGTH_FILTER.out.converted_bed)

    // Extract introns from merged transcripts
    INTRON_CHAIN_EXTRACTOR_MODELS(merged_transcripts)
    
    // Create full length contains fileron_chains_reads.join(INTRON_CHAIN_EXTRACTOR_MODELS.out.intron_chains_models, by: [0,1,2]).view()
    FULL_LENGTH_CONTAINS(INTRON_CHAIN_EXTRACTOR_READS.out.intron_chains_reads.join(INTRON_CHAIN_EXTRACTOR_MODELS.out.intron_chains_models, by: [0,1,2]))

    // Output full length unfiltered contains file
    COMBINE_UNFILTERED_FULL_LENGTH_CONTAINS(FULL_LENGTH_CONTAINS.out.full_length_contains.map { _chromosome_id, genome_id, contains, _stats, chromosome_count -> tuple(groupKey(genome_id, chromosome_count.toInteger()), contains) }.groupTuple())

    // Output full length stats
    COMBINE_FULL_LENGTH_STATS(FULL_LENGTH_CONTAINS.out.full_length_contains.map { _chromosome_id, genome_id, _full_length_contains, full_length_contains_stats, chromosome_count -> tuple(groupKey(genome_id, chromosome_count.toInteger())), full_length_contains_stats }.groupTuple())

    // Filter gtf for at least the specified number of full length reads
    FILTER_GTF_FOR_FULL_LENGTH(merged_transcripts.map { genome_id, gtf, _contains, chromosome_count, chromosome_id -> tuple(chromosome_id, genome_id, chromosome_count, gtf) }.join(FULL_LENGTH_CONTAINS.out.full_length_contains.map  { chromosome_id, genome_id, full_length_contains, _full_length_contains_stats, chromosome_count -> tuple(chromosome_id, genome_id, chromosome_count, full_length_contains) }, by : [0,1,2]))
    
    // Filter reads for full length reads
    FILTER_BAM_FOR_FULL_LENGTH(merged_reads.map { bam, _bai, genome_id, chromosome_id, chromosome_count -> tuple(chromosome_id, genome_id, chromosome_count, bam) }.join(FULL_LENGTH_CONTAINS.out.full_length_contains.map  { chromosome_id, genome_id, full_length_contains, _full_length_contains_stats, chromosome_count -> tuple(chromosome_id, genome_id, chromosome_count, full_length_contains) }, by : [0,1,2]))
    
    // Combine filtered GTF
    COMBINE_FILTERED_MODELS(FILTER_GTF_FOR_FULL_LENGTH.out.full_length_filtered_gtf.map { genome_id, chromosome_count, gtf, contains -> tuple(groupKey(genome_id, chromosome_count.toInteger()), gtf, contains)}.groupTuple())
    
    // Combine full length reads BAM
    COMBINE_FULL_LENGTH_READS(FILTER_BAM_FOR_FULL_LENGTH.out.full_length_filtered_bam.map { genome_id, chromosome_count, bam, contains -> tuple(groupKey(genome_id, chromosome_count.toInteger()), bam, contains)}.groupTuple())

    emit:
    FULL_LENGTH_CONTAINS.out.full_length_contains
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
