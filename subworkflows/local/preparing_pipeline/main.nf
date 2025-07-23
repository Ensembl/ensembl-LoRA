
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
        fasta: it[1] == "fasta" ? [it, 'fasta'] : null
        fastq: it[1] == "fastq" ? [it, 'fastq'] : null
        bam: it[1] == "bam" ? [it, 'bam'] : null
    }
    
    emit:
    fasta_ch = split_ch.fasta
    fastq_ch = split_ch.fastq
    bam_ch = split_ch.bam.map { it -> [it[0], it[3], it[4], it[5]] }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
