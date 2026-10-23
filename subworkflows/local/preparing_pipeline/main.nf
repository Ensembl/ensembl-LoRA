
//
// Subworkflow for the preparation of the LoRA pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROCESS_SAMPLE_SHEETS } from '../../../modules/local/process_sample_sheets/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARING_PIPELINE {
    take:
    read_sample_sheet
    reference_sample_sheet

    main:
    // Process the input sample sheets
    PROCESS_SAMPLE_SHEETS(read_sample_sheet, reference_sample_sheet)
    
    // Branch the input channel based on file type
    split_ch = PROCESS_SAMPLE_SHEETS.out.processed_sample_sheets.splitCsv().branch { 
        fasta: it[1] == "fasta"
        fastq: it[1] == "fastq"
        bam: it[1] == "bam"
    }
    
    emit:
    fasta_ch = split_ch.fasta.map { it -> [it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[9], it[0], it[7], it[8]] }
    fastq_ch = split_ch.fastq.map { it -> [it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[9], it[0], it[7], it[8]] }
    bam_ch = split_ch.bam.map { it -> [it[0], it[3], it[4], it[5]] }
    sqanti3_reference_ch = PROCESS_SAMPLE_SHEETS.out.sqanti3_reference_sheets.splitCsv()
    sample_metadata = PROCESS_SAMPLE_SHEETS.out.sample_metadata
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
