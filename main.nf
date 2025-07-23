#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LoRA: Long-read RNA-seq Analysis pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : <TBA>
    Author: Ryan Merritt
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARING_PIPELINE       } from './subworkflows/local/preparing_pipeline/main.nf'
include { PREPARE_FAST_FILES       } from './subworkflows/local/prepare_fast_files/main.nf'
include { MAPPING_READS            } from './subworkflows/local/mapping/main.nf'
include { SAMTOOLS                 } from './subworkflows/local/samtools/main.nf'
include { BAM_TO_GTF               } from './subworkflows/local/bam_to_gtf/main.nf'
include { PROVIDED_GTF_PREPARATION } from './subworkflows/local/provided_gtf_preparation/main.nf'
include { MERGING_TRANSCRIPTS      } from './subworkflows/local/merging_transcripts/main.nf'
include { FULL_LENGTH_FILTER } from './subworkflows/local/full_length_filter/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Run Preparing_pipeline SubWorkflow
    PREPARING_PIPELINE("${params.sample_sheet}", "${params.reference_sheet}")

    // Run Preparation SubWorkflow
    PREPARE_FAST_FILES(PREPARING_PIPELINE.out.fastq_ch, PREPARING_PIPELINE.out.fasta_ch)

    // Run Mapping SubWorkflow
    MAPPING_READS(PREPARE_FAST_FILES.out.fasta_files)

    // Run Samtools SubWorkflow
    SAMTOOLS(MAPPING_READS.out.mapped_reads, PREPARING_PIPELINE.out.bam_ch)

    // Run Bam to GTF SubWorkflow
    BAM_TO_GTF(SAMTOOLS.out.prepared_bam, SAMTOOLS.out.chromosome_bam)

    // Run Merging Transcripts SubWorkflow
    MERGING_TRANSCRIPTS(BAM_TO_GTF.out.prepared_gtf)

    // Filtering reads based on full length reads
    FULL_LENGTH_FILTER(BAM_TO_GTF.out.chromosome_bam, MERGING_TRANSCRIPTS.out.transcript_models_per_chromosome)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/