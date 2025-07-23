//
// Subworkflow for the preparation of fasta and fastq files for the LoRA pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AVERAGE_FASTQ_READ_LENGTH } from "../../../modules/local/average_fastq_read_length/main.nf"
include { AVERAGE_FASTA_READ_LENGTH } from "../../../modules/local/average_fasta_read_length/main.nf"
include { EDIT_FASTQ_IDS            } from "../../../modules/local/edit_fastq_read_ids/main.nf"
include { EDIT_FASTA_IDS            } from '../../../modules/local/edit_fasta_read_ids/main.nf'
include { FASTQ_TO_FASTA            } from '../../../modules/local/fastq_to_fasta/main.nf'

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
    // Get the average read length
    AVERAGE_FASTQ_READ_LENGTH(fastq_files)
    AVERAGE_FASTA_READ_LENGTH(fasta_files)

    // Edit read ids
    EDIT_FASTQ_IDS(AVERAGE_FASTQ_READ_LENGTH.out.passed)
    EDIT_FASTA_IDS(AVERAGE_FASTA_READ_LENGTH.out.passed)

    // Convert FASTQ to FASTA
    FASTQ_TO_FASTA(EDIT_FASTQ_IDS.out.edited_fastq)

    emit:
    fasta_files = FASTQ_TO_FASTA.out.converted_fasta.mix(EDIT_FASTA_IDS.out.edited_fasta)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
