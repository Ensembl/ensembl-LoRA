#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LoRA: Long-read RNA-seq Analysis pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : <TBA>
    Author: Ryan Merritt
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { LORA_MAIN } from './workflows/lora.nf'
include { LORA_INITIALISATION } from './subworkflows/local/lora_initialisation/main.nf'

workflow {
    LORA_INITIALISATION(
        workflow,
        params.validate_params,
        params.monochrome_logs,
        params.help,
        params.help_full,
        params.show_hidden,
        args
    )

    sample_sheet = params.sample_sheet ?: params.input
    reference_sheet = params.reference_sheet

    LORA_MAIN(sample_sheet, reference_sheet)
}
