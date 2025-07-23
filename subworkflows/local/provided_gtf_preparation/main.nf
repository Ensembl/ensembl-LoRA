//
// Subworkflow for the preparation of provided GTFs for the LoRA pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MERGE_PROVIDED_GTF           } from '../../../modules/local/merge_provided_gtf/main.nf'
include { SPLIT_GTF_INTO_BINS          } from '../../../modules/local/split_gtf_into_bins/main.nf'
include { SPLIT_PROVIDED_GTF_CHANNEL   } from '../../../modules/local/split_provided_gtf_channel/main.nf'
include { SORT_PROVIDED_GTF            } from '../../../modules/local/sort_provided_gtf/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROVIDED_GTF_PREPARATION {
    take:
    provided_gtf

    main:
    // Merge provided gtfs if provided
    MERGE_PROVIDED_GTF(provided_gtf.map { genome_id, gtf, individual_bin_count, total_bin_count, gtf_count -> tuple(groupKey(genome_id, gtf_count.toInteger()), gtf, individual_bin_count, total_bin_count) }.groupTuple().map { genome_id, gtf, individual_bin_count, total_bin_count -> tuple(gtf[0], genome_id, individual_bin_count[0], total_bin_count[0]) })

    // Split the provided chromosome gtf files into bins
    SPLIT_GTF_INTO_BINS(MERGE_PROVIDED_GTF.out.merged_fixed_gtf)

    // Split SPLIT_GTF_INTO_BINS output into separate channels
    SPLIT_PROVIDED_GTF_CHANNEL(SPLIT_GTF_INTO_BINS.out.binned_gtf)

    // sort the split gtf files
    SORT_PROVIDED_GTF(SPLIT_PROVIDED_GTF_CHANNEL.out.prepared_gtf_sample_sheet.splitCsv())

    emit:
    provided_gtfs = SORT_PROVIDED_GTF.out.sorted_gtf
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/