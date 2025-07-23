//
// Subworkflow for the alignment of fasta files for the LoRA pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PB_MINIMAP2 } from '../../../modules/local/pb_minimap2/main.nf'
include { NP_MINIMAP2 } from '../../../modules/local/np_minimap2/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAPPING_READS {
    take:
    fasta_files

    main:
    // Branch the fasta files based on sequencing source
    split_ch = fasta_files.branch {
        pacbio: it[1] == "pacbio"
        nanopore: it[1] == "nanopore"
    }

    // Align pacbio reads
    PB_MINIMAP2(split_ch.pacbio)

    // Align nanopore reads
    NP_MINIMAP2(split_ch.nanopore)

    emit:
    mapped_reads = PB_MINIMAP2.out.mapped_reads_pb.mix(NP_MINIMAP2.out.mapped_reads_np)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
