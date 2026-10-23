//
// Optional SQANTI3 reporting branch for LoRA transcript models
//

include { PREPARE_SQANTI3_REFERENCE } from '../../../modules/local/prepare_sqanti3_reference/main.nf'
include { PREPARE_SQANTI3_CHROMOSOME_REFERENCE } from '../../../modules/local/prepare_sqanti3_chromosome_reference/main.nf'
include { WRITE_SQANTI3_INPUT_MANIFEST } from '../../../modules/local/write_sqanti3_input_manifest/main.nf'
include { VALIDATE_SQANTI3_INPUT_MANIFESTS } from '../../../modules/local/validate_sqanti3_input_manifests/main.nf'
include { SQANTI3_QC                } from '../../../modules/local/sqanti3_qc/main.nf'
include { MERGE_SQANTI3_OUTPUTS     } from '../../../modules/local/merge_sqanti3_outputs/main.nf'

workflow SQANTI3 {
    take:
    sqanti3_references
    prefiltered_models
    filtered_models
    prefiltered_models_per_chromosome
    filtered_models_per_chromosome

    main:
    PREPARE_SQANTI3_REFERENCE(sqanti3_references)

    sqanti3_by_chromosome = params.sqanti3_by_chromosome.toString().toLowerCase() in ['true', '1', 'yes']

    sqanti3_stage_models = prefiltered_models
        .map { genome_id, isoforms -> tuple("prefiltered", genome_id, isoforms) }
        .mix(
            filtered_models
                .map { genome_id, isoforms, _contains -> tuple("postfiltered", genome_id, isoforms) }
        )

    if (sqanti3_by_chromosome) {
        sqanti3_chromosome_models = prefiltered_models_per_chromosome
            .map { genome_id, gtf, _contains, _chromosome_count, chromosome_id ->
                tuple(genome_id, "prefiltered", chromosome_id, gtf)
            }
            .mix(
                filtered_models_per_chromosome
                    .map { genome_id, _chromosome_count, chromosome_id, gtf, _contains ->
                        tuple(genome_id, "postfiltered", chromosome_id, gtf)
                    }
            )
            .filter { genome_id, stage, chromosome_id, isoforms ->
                def file = new File(isoforms.toString())
                file.exists() && file.length() > 0
            }

        sqanti3_chromosome_reference_inputs = sqanti3_chromosome_models
            .combine(PREPARE_SQANTI3_REFERENCE.out.prepared_reference, by: 0)
            .map { genome_id, stage, chromosome_id, isoforms, reference_fasta, reference_gtf, polya_peak, polya_source, cage_peak, cage_source ->
                tuple(genome_id, stage, chromosome_id, isoforms, reference_fasta, reference_gtf, polya_peak, polya_source, cage_peak, cage_source)
            }

        PREPARE_SQANTI3_CHROMOSOME_REFERENCE(sqanti3_chromosome_reference_inputs)

        sqanti3_candidate_inputs = PREPARE_SQANTI3_CHROMOSOME_REFERENCE.out.chromosome_reference
            .map { stage, genome_id, chromosome_id, isoforms, reference_fasta, reference_gtf, polya_peak, polya_source, cage_peak, cage_source, _chromosome_map ->
                tuple(stage, genome_id, chromosome_id, isoforms, reference_fasta, reference_gtf, polya_peak, polya_source, cage_peak, cage_source)
            }
    } else {
        sqanti3_chunk_models = sqanti3_stage_models
            .map { stage, genome_id, isoforms -> tuple(genome_id, stage, "full", isoforms) }

        sqanti3_candidate_inputs = sqanti3_chunk_models
            .combine(PREPARE_SQANTI3_REFERENCE.out.prepared_reference, by: 0)
            .map { genome_id, stage, chunk_id, isoforms, reference_fasta, reference_gtf, polya_peak, polya_source, cage_peak, cage_source ->
                tuple(stage, genome_id, chunk_id, isoforms, reference_fasta, reference_gtf, polya_peak, polya_source, cage_peak, cage_source)
            }
    }

    WRITE_SQANTI3_INPUT_MANIFEST(sqanti3_candidate_inputs)

    VALIDATE_SQANTI3_INPUT_MANIFESTS(WRITE_SQANTI3_INPUT_MANIFEST.out.sqanti3_input_manifest.collect())

    sqanti3_qc_inputs = VALIDATE_SQANTI3_INPUT_MANIFESTS.out.sqanti3_qc_inputs
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.stage,
                row.genome_id,
                row.chunk_id,
                file(row.isoforms),
                file(row.reference_fasta),
                file(row.reference_gtf),
                file(row.polya_peak),
                row.polya_source,
                file(row.cage_peak),
                row.cage_source
            )
        }

    SQANTI3_QC(sqanti3_qc_inputs)

    sqanti3_merge_inputs = SQANTI3_QC.out.sqanti3_results
        .groupTuple(by: [0, 1])
        .map { stage, genome_id, chunk_ids, chunk_dirs ->
            tuple(stage, genome_id, chunk_ids, chunk_dirs)
        }

    MERGE_SQANTI3_OUTPUTS(sqanti3_merge_inputs)

    emit:
    sqanti3_results = MERGE_SQANTI3_OUTPUTS.out.sqanti3_results
}
