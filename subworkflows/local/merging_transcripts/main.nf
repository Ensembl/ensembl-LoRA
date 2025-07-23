//
// Subworkflow for the merging of transcripts for the LoRA pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FIRST_TMERGE                    } from '../../../modules/local/first_tmerge/main.nf'
include { FIX_FIRST_TMERGE                } from '../../../modules/local/fix_first_tmerge/main.nf'
include { FIRST_TMERGE_TO_GTF             } from '../../../modules/local/first_tmerge_to_gtf/main.nf'
include { COMBINE_FIRST_TMERGE            } from '../../../modules/local/combine_first_tmerge/main.nf'
include { EXTRACT_CHROMOSOME_IDS_FROM_GTF } from '../../../modules/local/extract_chromosome_ids_from_gtf/main.nf'
include { SPLIT_GTF_BY_CHROMOSOME         } from '../../../modules/local/split_gtf_by_chromosome/main.nf'
include { SECOND_TMERGE                   } from '../../../modules/local/second_tmerge/main.nf'
include { FIX_SECOND_TMERGE               } from '../../../modules/local/fix_second_tmerge/main.nf'
include { SECOND_TMERGE_TO_GTF            } from '../../../modules/local/second_tmerge_to_gtf/main.nf'
include { BUILD_CONTAINS_FILES            } from '../../../modules/local/build_contains_files/main.nf'
include { COMBINE_TRANSCRIPT_MODELS       } from '../../../modules/local/combine_transcript_models/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MERGING_TRANSCRIPTS {
    take:
    transcript_models

    main:
    // Run the first tmerge process
    FIRST_TMERGE(transcript_models)

    // Fix the first tmerge output
    FIX_FIRST_TMERGE(FIRST_TMERGE.out.first_tmerge)

    // Convert the first tmerge output to a normal gtf file
    FIRST_TMERGE_TO_GTF(FIX_FIRST_TMERGE.out.first_tmerge_fixed)

    // Combine gtfs by chromosome
    COMBINE_FIRST_TMERGE(FIRST_TMERGE_TO_GTF.out.first_tmerge_gtf.map { chromosome_id, gtf, contains, genome_id, chromosome_count, bin_count -> tuple(groupKey(chromosome_id, bin_count.toInteger()), gtf, contains, genome_id, chromosome_count)}.groupTuple().map { chromosome_id, gtf, contains, genome_id, chromosome_count -> [genome_id[0], gtf, contains, chromosome_count[0], chromosome_id] })
    
    // Running the final tmerge
    SECOND_TMERGE(COMBINE_FIRST_TMERGE.out.combined_bins)

    // Fix the second tmerge output
    FIX_SECOND_TMERGE(SECOND_TMERGE.out.second_tmerge_output)

    // Convert the final tmerge output to a normal gtf file
    SECOND_TMERGE_TO_GTF(FIX_SECOND_TMERGE.out.second_tmerge_fixed)

    // Building the contains files
    BUILD_CONTAINS_FILES(SECOND_TMERGE_TO_GTF.out.retmerge_to_gtf_output)
    
    // Combine the transcript models and the contains files
    COMBINE_TRANSCRIPT_MODELS(BUILD_CONTAINS_FILES.out.combined_contains.map { genome_id, gtf, contains, chromosome_count, chromosome_id -> tuple(groupKey(genome_id, chromosome_count.toInteger()), gtf, contains, chromosome_id) }.groupTuple())

    emit:
    transcript_models_per_chromosome = BUILD_CONTAINS_FILES.out.combined_contains
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
