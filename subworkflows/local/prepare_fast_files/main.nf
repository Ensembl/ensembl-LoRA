//
// Subworkflow for the preparation of fasta and fastq files for the LoRA pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EDIT_FASTQ_IDS            } from "../../../modules/local/edit_fastq_read_ids/main.nf"
include { EDIT_FASTA_IDS            } from '../../../modules/local/edit_fasta_read_ids/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_FAST_FILES {
    take:
    fastq_files
    fasta_files

    main:
    // Validate read lengths while standardising read ids
    EDIT_FASTQ_IDS(fastq_files)
    EDIT_FASTA_IDS(fasta_files)

    emit:
    fasta_files = EDIT_FASTQ_IDS.out.edited_fasta.mix(EDIT_FASTA_IDS.out.edited_fasta)
    read_id_maps = EDIT_FASTQ_IDS.out.read_id_map.mix(EDIT_FASTA_IDS.out.read_id_map)
    read_summaries = EDIT_FASTQ_IDS.out.read_summary.mix(EDIT_FASTA_IDS.out.read_summary)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
