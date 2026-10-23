/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LoRA workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def checkpointStages() {
    return ['read_prep', 'mapped_bam', 'bam_to_gtf', 'merged_models', 'full_length', 'sqanti3', 'dashboard']
}

def implementedCheckpointStages() {
    return ['read_prep', 'mapped_bam', 'merged_models', 'full_length', 'sqanti3', 'dashboard']
}

def asBool(value) {
    return value.toString().toLowerCase() in ['true', '1', 'yes']
}

def resolveOutdirPath(outdir, launchDir) {
    def raw = outdir.toString()
    def file = new File(raw)
    if (file.isAbsolute()) {
        return file.canonicalPath
    }
    return new File(launchDir.toString(), raw).canonicalPath
}

def nonEmptyFile(path) {
    def file = new File(path.toString())
    return file.exists() && file.isFile() && file.length() > 0
}

def genomeDirs(outdir) {
    def root = new File(outdir.toString())
    def excluded = ['dashboard', 'pipeline_info', 'summary', 'work', 'tmp', 'checkpoints'] as Set
    if (!root.exists()) {
        return []
    }
    return root.listFiles()
        .findAll { it.isDirectory() && !excluded.contains(it.name) && !it.name.startsWith('.') }
        .collect { it.name }
        .sort()
}

def checkpointOutdirMessage(outdir) {
    def root = new File(outdir.toString())
    def lines = []
    lines << "Scanned checkpoint outdir: ${root.path}"
    lines << "Expected genome directories like: ${root.path}/<genome_id>/"
    if (!root.exists()) {
        lines << "The checkpoint outdir does not exist."
    } else if (!root.isDirectory()) {
        lines << "The checkpoint outdir exists but is not a directory."
    } else {
        def entries = root.listFiles()?.collect { it.name }?.sort()?.take(30) ?: []
        lines << "Found top-level entries: ${entries ? entries.join(', ') : '(none)'}"
    }
    if (root.path.startsWith('/nfs/')) {
        lines << "This looks like an EBI /nfs path. If the run used /gpfs, try the matching /gpfs path."
    } else if (root.path.startsWith('/gpfs/')) {
        lines << "This looks like an EBI /gpfs path. If the run used /nfs, check whether the same outputs exist through that mount."
    }
    lines << "Run: bin/lora_checkpoint_status.py --outdir ${root.path} --run_id ${params.runID}"
    return lines.join('\n')
}

def checkpointFile(outdir, runID, genome_id, stage, kind) {
    def exact = [
        read_summary: "${outdir}/${genome_id}/reads/read_summaries/${runID}_${genome_id}_read_summary.tsv",
        prefiltered_gtf: "${outdir}/${genome_id}/all_models/merged/${runID}_${genome_id}.gtf",
        prefiltered_contains: "${outdir}/${genome_id}/all_models/merged/${runID}_${genome_id}_contains.tsv",
        postfiltered_gtf: "${outdir}/${genome_id}/full_length/models/merged/${runID}_${genome_id}_filtered_transcripts.gtf",
        postfiltered_contains: "${outdir}/${genome_id}/full_length/models/merged/${runID}_${genome_id}_filtered_contains.tsv",
        full_length_stats: "${outdir}/${genome_id}/full_length/stats/${runID}_${genome_id}_full_length_stats.tsv",
        sample_metadata: "${outdir}/summary/sample_metadata.csv",
        sqanti3_dir: "${outdir}/${genome_id}/sqanti3",
        dashboard_data: "${outdir}/${genome_id}/summary/${runID}_${genome_id}_dashboard_data.json"
    ]
    return exact[kind]
}

def preparedFastaPublished(outdir, genome_id) {
    def genomeDir = new File("${outdir}/${genome_id}")
    def dirs = [
        'reads/prepared_fasta',
        'reads/prepared',
        'reads/fasta',
        'reads/renamed',
        'reads/all_reads/fasta',
        'reads/all_reads/renamed'
    ]
    return dirs.any { rel ->
        def dir = new File(genomeDir, rel)
        dir.exists() && dir.listFiles()?.any { it.isFile() && it.name.endsWith('.fasta.gz') && it.length() > 0 }
    }
}

def checkpointStatus(outdir, runID, genome_id, stage) {
    if (stage == 'read_prep') {
        if (nonEmptyFile(checkpointFile(outdir, runID, genome_id, stage, 'read_summary')) && preparedFastaPublished(outdir, genome_id)) {
            return 'complete'
        }
        return nonEmptyFile(checkpointFile(outdir, runID, genome_id, stage, 'read_summary')) ? 'not_resumable' : 'missing'
    }
    if (stage == 'mapped_bam') {
        def bam = "${outdir}/${genome_id}/reads/all_reads/merged/${runID}_${genome_id}.bam"
        def bai = "${outdir}/${genome_id}/reads/all_reads/merged/${runID}_${genome_id}.bam.csi"
        return (nonEmptyFile(bam) && nonEmptyFile(bai)) ? 'complete' : 'missing'
    }
    if (stage == 'bam_to_gtf') {
        def canonicalDir = new File("${outdir}/${genome_id}/reads/canonical/chromosomes")
        def modelDir = new File("${outdir}/${genome_id}/all_models/chromosomes")
        def hasCanonical = canonicalDir.exists() && canonicalDir.listFiles()?.any { it.isFile() && it.name.endsWith('.canonical.bam') && it.length() > 0 }
        def hasChromModels = modelDir.exists() && modelDir.listFiles()?.any { it.isFile() && it.name.endsWith('.gtf') && it.length() > 0 }
        return (hasCanonical || hasChromModels) ? 'not_resumable' : 'missing'
    }
    if (stage == 'merged_models') {
        def chromDir = new File("${outdir}/${genome_id}/all_models/chromosomes")
        def canonicalDir = new File("${outdir}/${genome_id}/reads/canonical/chromosomes")
        def hasChromModels = chromDir.exists() && chromDir.listFiles()?.any { it.isFile() && it.name.endsWith('.combined.gtf') && it.length() > 0 }
        def hasCanonical = canonicalDir.exists() && canonicalDir.listFiles()?.any { it.isFile() && it.name.endsWith('.canonical.bam') && it.length() > 0 }
        return (
            nonEmptyFile(checkpointFile(outdir, runID, genome_id, stage, 'prefiltered_gtf')) &&
            nonEmptyFile(checkpointFile(outdir, runID, genome_id, stage, 'prefiltered_contains')) &&
            hasChromModels &&
            hasCanonical
        ) ? 'complete' : 'missing'
    }
    if (stage == 'full_length') {
        return (
            nonEmptyFile(checkpointFile(outdir, runID, genome_id, stage, 'prefiltered_gtf')) &&
            nonEmptyFile(checkpointFile(outdir, runID, genome_id, stage, 'prefiltered_contains')) &&
            nonEmptyFile(checkpointFile(outdir, runID, genome_id, stage, 'postfiltered_gtf')) &&
            nonEmptyFile(checkpointFile(outdir, runID, genome_id, stage, 'postfiltered_contains')) &&
            nonEmptyFile(checkpointFile(outdir, runID, genome_id, stage, 'full_length_stats')) &&
            nonEmptyFile(checkpointFile(outdir, runID, genome_id, stage, 'read_summary')) &&
            nonEmptyFile(checkpointFile(outdir, runID, genome_id, stage, 'sample_metadata'))
        ) ? 'complete' : 'missing'
    }
    if (stage == 'sqanti3') {
        def root = new File(checkpointFile(outdir, runID, genome_id, stage, 'sqanti3_dir'))
        if (!root.exists()) {
            return 'skipped'
        }
        def pre = new File("${root}/prefiltered/${runID}_${genome_id}_prefiltered_classification.txt")
        def post = new File("${root}/postfiltered/${runID}_${genome_id}_postfiltered_classification.txt")
        return (pre.exists() && pre.length() > 0 && post.exists() && post.length() > 0) ? 'complete' : 'missing'
    }
    if (stage == 'dashboard') {
        def genomeDashboard = new File("${outdir}/${genome_id}/dashboard/dashboard_data.json")
        def centralDashboard = new File("${outdir}/dashboard/${genome_id}/dashboard_data.json")
        return ((genomeDashboard.exists() && genomeDashboard.length() > 0) || (centralDashboard.exists() && centralDashboard.length() > 0)) ? 'complete' : 'missing'
    }
    return 'missing'
}

def latestCommonCheckpoint(outdir, runID, genomes) {
    def latest = 'none'
    checkpointStages().each { stage ->
        if (implementedCheckpointStages().contains(stage) && genomes && genomes.every { checkpointStatus(outdir, runID, it, stage) == 'complete' }) {
            latest = stage
        }
    }
    return latest
}

def requiredCheckpointsForResume(stage) {
    if (stage == 'sqanti3') {
        return ['full_length', 'sqanti3']
    }
    if (stage == 'dashboard') {
        return ['full_length', 'dashboard']
    }
    return [stage]
}

def resolveCheckpointPlan(params, checkpoint_outdir) {
    def requested = params.resume_from.toString()
    def scope = params.resume_scope.toString()
    if (requested == 'none') {
        return [:]
    }
    if (requested == 'bam_to_gtf') {
        throw new IllegalArgumentException("bam_to_gtf is audit-only in this LoRA version. Resume from mapped_bam to rerun BAM-to-GTF, or from merged_models if model merging completed.")
    }
    def genomes = genomeDirs(checkpoint_outdir)
    if (!genomes) {
        throw new IllegalArgumentException("No genome output directories found; cannot use published-output checkpoint resume.\n${checkpointOutdirMessage(checkpoint_outdir)}")
    }
    if (scope == 'common') {
        def resolved = requested == 'auto' ? latestCommonCheckpoint(checkpoint_outdir, params.runID, genomes) : requested
        validateResolvedCheckpoint(params, checkpoint_outdir, resolved, genomes)
        return genomes.collectEntries { genome_id -> [(genome_id): resolved] }
    }
    if (scope == 'per_genome') {
        return genomes.collectEntries { genome_id ->
            def resolved = requested == 'auto' ? latestGenomeCheckpoint(checkpoint_outdir, params.runID, genome_id) : requested
            validateResolvedCheckpoint(params, checkpoint_outdir, resolved, [genome_id])
            [(genome_id): resolved]
        }
    }
    throw new IllegalArgumentException("Unsupported resume_scope='${scope}'. Use common or per_genome.")
}

def latestGenomeCheckpoint(outdir, runID, genome_id) {
    def latest = 'none'
    checkpointStages().each { stage ->
        if (implementedCheckpointStages().contains(stage) && checkpointStatus(outdir, runID, genome_id, stage) == 'complete') {
            latest = stage
        }
    }
    return latest
}

def validateResolvedCheckpoint(params, checkpoint_outdir, resolved, genomes) {
    if (resolved == 'bam_to_gtf') {
        throw new IllegalArgumentException("bam_to_gtf is audit-only in this LoRA version. Resume from mapped_bam to rerun BAM-to-GTF, or from merged_models if model merging completed.")
    }
    if (!implementedCheckpointStages().contains(resolved)) {
        throw new IllegalArgumentException("Checkpoint resume_from='${resolved}' is not resumable from published outputs. Implemented stages: ${implementedCheckpointStages().join(', ')}. Run bin/lora_checkpoint_status.py --outdir ${checkpoint_outdir} --run_id ${params.runID}.")
    }
    requiredCheckpointsForResume(resolved).each { required_stage ->
        def incomplete = genomes.findAll { checkpointStatus(checkpoint_outdir, params.runID, it, required_stage) != 'complete' }
        if (incomplete) {
            throw new IllegalArgumentException("Checkpoint '${resolved}' requires '${required_stage}', but it is incomplete for genome(s): ${incomplete.join(', ')}. Run bin/lora_checkpoint_status.py --outdir ${checkpoint_outdir} --run_id ${params.runID}.")
        }
    }
    return resolved
}

include { PREPARING_PIPELINE       } from '../subworkflows/local/preparing_pipeline/main.nf'
include { PREPARE_FAST_FILES       } from '../subworkflows/local/prepare_fast_files/main.nf'
include { MAPPING_READS            } from '../subworkflows/local/mapping/main.nf'
include { MAPPING_READS as MAPPING_READS_CHECKPOINT } from '../subworkflows/local/mapping/main.nf'
include { SAMTOOLS                 } from '../subworkflows/local/samtools/main.nf'
include { SAMTOOLS as SAMTOOLS_CHECKPOINT } from '../subworkflows/local/samtools/main.nf'
include { SAMTOOLS_FROM_MERGED_BAM } from '../subworkflows/local/samtools_from_merged_bam/main.nf'
include { BAM_TO_GTF               } from '../subworkflows/local/bam_to_gtf/main.nf'
include { BAM_TO_GTF as BAM_TO_GTF_CHECKPOINT } from '../subworkflows/local/bam_to_gtf/main.nf'
include { PROVIDED_GTF_PREPARATION } from '../subworkflows/local/provided_gtf_preparation/main.nf'
include { MERGING_TRANSCRIPTS      } from '../subworkflows/local/merging_transcripts/main.nf'
include { MERGING_TRANSCRIPTS as MERGING_TRANSCRIPTS_CHECKPOINT } from '../subworkflows/local/merging_transcripts/main.nf'
include { FULL_LENGTH_FILTER       } from '../subworkflows/local/full_length_filter/main.nf'
include { FULL_LENGTH_FILTER as FULL_LENGTH_FILTER_CHECKPOINT } from '../subworkflows/local/full_length_filter/main.nf'
include { SQANTI3                  } from '../subworkflows/local/sqanti3/main.nf'
include { BUILD_LORA_DASHBOARD     } from '../modules/local/build_lora_dashboard/main.nf'
include { BUILD_LORA_DASHBOARD_INDEX } from '../modules/local/build_lora_dashboard_index/main.nf'
include { COMBINE_READ_ID_MAPS     } from '../modules/local/combine_read_id_maps/main.nf'
include { COMBINE_READ_SUMMARIES   } from '../modules/local/combine_read_summaries/main.nf'
include { PREPARE_SQANTI3_DASHBOARD_INPUT } from '../modules/local/prepare_sqanti3_dashboard_input/main.nf'
include { PREPARE_EMPTY_SQANTI3_DASHBOARD_INPUT } from '../modules/local/prepare_empty_sqanti3_dashboard_input/main.nf'
include { COUNT_LORA_SAMPLES       } from '../modules/local/count_lora_samples/main.nf'
include { COUNT_LORA_INPUT_READS   } from '../modules/local/count_lora_input_reads/main.nf'
include { COUNT_LORA_MODELS        } from '../modules/local/count_lora_models/main.nf'
include { COUNT_LORA_SQANTI3       } from '../modules/local/count_lora_sqanti3/main.nf'
include { CHECK_LORA_BAM_STATS     } from '../modules/local/check_lora_bam_stats/main.nf'
include { BUILD_LORA_DATA_CHECKS   } from '../modules/local/build_lora_data_checks/main.nf'
include { BUILD_LORA_DASHBOARD_MANIFESTS } from '../modules/local/build_lora_dashboard_manifests/main.nf'
include { SUMMARISE_LORA_DASHBOARD } from '../modules/local/summarise_lora_dashboard/main.nf'
include { LOAD_LORA_CHECKPOINTS    } from '../modules/local/load_lora_checkpoints/main.nf'

workflow LORA {

    take:
    sample_sheet
    reference_sheet

    main:
    run_sqanti3 = asBool(params.run_sqanti3)
    run_dashboard = asBool(params.dashboard)
    combine_read_id_maps = asBool(params.combine_read_id_maps)
    checkpoint_outdir = resolveOutdirPath(params.outdir, workflow.launchDir)
    checkpoint_plan = resolveCheckpointPlan(params, checkpoint_outdir)
    checkpoint_resume = checkpoint_plan ? 'checkpoint' : 'none'
    checkpoint_genomes = checkpoint_plan.keySet().toList().sort()
    read_prep_genomes = checkpoint_plan.findAll { genome_id, stage -> stage == 'read_prep' }.keySet() as Set
    sqanti_run_genomes = checkpoint_plan.findAll { genome_id, stage -> stage in ['read_prep', 'mapped_bam', 'merged_models', 'full_length'] }.keySet() as Set
    sqanti_existing_genomes = checkpoint_plan.findAll { genome_id, stage -> stage in ['sqanti3', 'dashboard'] }.keySet() as Set
    dashboard_metric_checkpoint_genomes = checkpoint_plan.findAll { genome_id, stage -> stage == 'dashboard' }.keySet() as Set
    dashboard_only_resume = checkpoint_plan && checkpoint_plan.values().every { stage -> stage == 'dashboard' }

    PREPARING_PIPELINE(sample_sheet, reference_sheet)

    if (checkpoint_resume == 'none') {
        PREPARE_FAST_FILES(PREPARING_PIPELINE.out.fastq_ch, PREPARING_PIPELINE.out.fasta_ch)

        if (combine_read_id_maps) {
            COMBINE_READ_ID_MAPS(
                PREPARE_FAST_FILES.out.read_id_maps
                    .map { read_id_map, _seq_source, _reference, genome_id, raw_file_count, _species ->
                        tuple(groupKey(genome_id, raw_file_count.toInteger()), read_id_map)
                    }
                    .groupTuple()
                    .map { genome_id, read_id_maps -> tuple(genome_id, read_id_maps) }
            )
            read_id_maps_ch = COMBINE_READ_ID_MAPS.out.combined_read_id_map
        } else {
            read_id_maps_ch = Channel.empty()
        }

        COMBINE_READ_SUMMARIES(
            PREPARE_FAST_FILES.out.read_summaries
                .map { read_summary, _seq_source, _reference, genome_id, raw_file_count, _species ->
                    tuple(groupKey(genome_id, raw_file_count.toInteger()), read_summary)
                }
                .groupTuple()
                .map { genome_id, read_summaries -> tuple(genome_id, read_summaries) }
        )

        MAPPING_READS(PREPARE_FAST_FILES.out.fasta_files)

        SAMTOOLS(MAPPING_READS.out.mapped_reads, PREPARING_PIPELINE.out.bam_ch)

        BAM_TO_GTF(SAMTOOLS.out.prepared_bam, SAMTOOLS.out.chromosome_bam)

        MERGING_TRANSCRIPTS(BAM_TO_GTF.out.prepared_gtf)

        FULL_LENGTH_FILTER(BAM_TO_GTF.out.chromosome_bam, MERGING_TRANSCRIPTS.out.transcript_models_per_chromosome)

        prefiltered_models_ch = MERGING_TRANSCRIPTS.out.prefiltered_models
        prefiltered_contains_ch = MERGING_TRANSCRIPTS.out.prefiltered_contains
        prefiltered_models_per_chromosome_ch = MERGING_TRANSCRIPTS.out.transcript_models_per_chromosome
        filtered_models_ch = FULL_LENGTH_FILTER.out.filtered_models
        filtered_models_per_chromosome_ch = FULL_LENGTH_FILTER.out.filtered_models_per_chromosome
        full_length_stats_ch = FULL_LENGTH_FILTER.out.full_length_stats
        read_summaries_ch = COMBINE_READ_SUMMARIES.out.combined_read_summary
        dashboard_bam_stats_inputs_ch = SAMTOOLS.out.sorted_bam_with_index
            .map { genome_id, bam, bai, _reference, _raw_file_count -> tuple('sample', genome_id, bam.simpleName.toString(), bam, bai) }
            .mix(SAMTOOLS.out.merged_bam
                .map { bam, bai, _reference, genome_id -> tuple('merged', genome_id, genome_id, bam, bai) })
            .mix(SAMTOOLS.out.chromosome_bam
                .map { bam, bai, _reference, genome_id, chromosome_id, _chromosome_count, _bin_count -> tuple('chromosome', genome_id, chromosome_id, bam, bai) })
            .mix(BAM_TO_GTF.out.chromosome_bam
                .map { bam, bai, genome_id, chromosome_id, _chromosome_count -> tuple('canonical', genome_id, chromosome_id, bam, bai) })
        dashboard_genomes_ch = prefiltered_models_ch.map { genome_id, _prefiltered_gtf -> genome_id }
        dashboard_empty_sqanti3_genomes_ch = prefiltered_models_ch.map { genome_id, _prefiltered_gtf -> genome_id }
        dashboard_manifest_genomes_ch = prefiltered_models_ch.map { genome_id, _prefiltered_gtf -> genome_id }
        dashboard_prefiltered_models_ch = prefiltered_models_ch
        dashboard_prefiltered_contains_ch = prefiltered_contains_ch
        dashboard_filtered_models_ch = filtered_models_ch
        dashboard_full_length_stats_ch = full_length_stats_ch
        dashboard_read_summaries_ch = read_summaries_ch
        checkpoint_sample_counts_dashboard_ch = Channel.empty()
        checkpoint_input_read_counts_dashboard_ch = Channel.empty()
        checkpoint_prefiltered_model_counts_dashboard_metric_ch = Channel.empty()
        checkpoint_postfiltered_model_counts_dashboard_metric_ch = Channel.empty()
        checkpoint_prefiltered_sqanti3_counts_dashboard_metric_ch = Channel.empty()
        checkpoint_postfiltered_sqanti3_counts_dashboard_metric_ch = Channel.empty()
        checkpoint_full_length_stats_dashboard_metric_ch = Channel.empty()
        checkpoint_data_checks_dashboard_ch = Channel.empty()
        checkpoint_dashboard_bam_stats_inputs_ch = Channel.empty()
        sqanti3_existing_dashboard_ch = Channel.empty()
    } else {
        log.info "LoRA published-output checkpoint resume: ${checkpoint_plan.collect { genome_id, stage -> genome_id + ':' + stage }.join(', ')}"
        LOAD_LORA_CHECKPOINTS(checkpoint_outdir, PREPARING_PIPELINE.out.sample_metadata, file(reference_sheet))

        dashboard_genomes_ch = Channel.fromList(checkpoint_genomes)
        dashboard_empty_sqanti3_genomes_ch = Channel.fromList(checkpoint_genomes)
        dashboard_manifest_genomes_ch = Channel.fromList(checkpoint_genomes)

        checkpoint_read_prep_fasta_ch = LOAD_LORA_CHECKPOINTS.out.read_prep_fasta
            .splitCsv()
            .map { row -> tuple(file(row[0]), row[1], file(row[2]), row[3], row[4], row[5]) }

        checkpoint_mapped_bam_ch = LOAD_LORA_CHECKPOINTS.out.mapped_bam
            .splitCsv()
            .map { row -> tuple(file(row[0]), file(row[1]), file(row[2]), row[3]) }

        checkpoint_canonical_bam_ch = LOAD_LORA_CHECKPOINTS.out.canonical_bam
            .splitCsv()
            .map { row -> tuple(file(row[0]), file(row[1]), row[2], row[3], row[4]) }

        checkpoint_chromosome_models_ch = LOAD_LORA_CHECKPOINTS.out.chromosome_models
            .splitCsv()
            .map { row -> tuple(row[0], file(row[1]), file(row[2]), row[3], row[4]) }

        checkpoint_prefiltered_models_ch = LOAD_LORA_CHECKPOINTS.out.prefiltered_models
            .splitCsv()
            .map { row -> tuple(row[0], file(row[1])) }

        checkpoint_prefiltered_contains_ch = LOAD_LORA_CHECKPOINTS.out.prefiltered_contains
            .splitCsv()
            .map { row -> tuple(row[0], file(row[1])) }

        checkpoint_filtered_models_ch = LOAD_LORA_CHECKPOINTS.out.filtered_models
            .splitCsv()
            .map { row -> tuple(row[0], file(row[1]), file(row[2])) }

        checkpoint_filtered_chromosome_models_ch = LOAD_LORA_CHECKPOINTS.out.filtered_chromosome_models
            .splitCsv()
            .map { row -> tuple(row[0], row[3], row[4], file(row[1]), file(row[2])) }

        checkpoint_full_length_stats_ch = LOAD_LORA_CHECKPOINTS.out.full_length_stats
            .splitCsv()
            .map { row -> tuple(row[0], file(row[1])) }

        read_summaries_ch = LOAD_LORA_CHECKPOINTS.out.read_summaries
            .splitCsv()
            .map { row -> tuple(row[0], file(row[1])) }

        sqanti3_existing_dashboard_ch = LOAD_LORA_CHECKPOINTS.out.sqanti3_dashboard
            .splitCsv()
            .map { row -> tuple(row[0], row[1].toString().toLowerCase() == 'true', file(row[2])) }

        checkpoint_prefiltered_models_dashboard_ch = LOAD_LORA_CHECKPOINTS.out.prefiltered_models
            .splitCsv()
            .map { row -> tuple(row[0], file(row[1])) }

        checkpoint_prefiltered_contains_dashboard_ch = LOAD_LORA_CHECKPOINTS.out.prefiltered_contains
            .splitCsv()
            .map { row -> tuple(row[0], file(row[1])) }

        checkpoint_filtered_models_dashboard_ch = LOAD_LORA_CHECKPOINTS.out.filtered_models
            .splitCsv()
            .map { row -> tuple(row[0], file(row[1]), file(row[2])) }

        checkpoint_full_length_stats_dashboard_ch = LOAD_LORA_CHECKPOINTS.out.full_length_stats
            .splitCsv()
            .map { row -> tuple(row[0], file(row[1])) }

        checkpoint_read_summaries_dashboard_ch = LOAD_LORA_CHECKPOINTS.out.read_summaries
            .splitCsv()
            .map { row -> tuple(row[0], file(row[1])) }

        checkpoint_sample_counts_dashboard_ch = LOAD_LORA_CHECKPOINTS.out.sample_counts
            .splitCsv()
            .filter { row -> dashboard_metric_checkpoint_genomes.contains(row[0]) }
            .map { row -> tuple(row[0], file(row[1])) }

        checkpoint_input_read_counts_dashboard_ch = LOAD_LORA_CHECKPOINTS.out.input_read_counts
            .splitCsv()
            .filter { row -> dashboard_metric_checkpoint_genomes.contains(row[0]) }
            .map { row -> tuple(row[0], file(row[1])) }

        checkpoint_prefiltered_model_counts_dashboard_metric_ch = LOAD_LORA_CHECKPOINTS.out.prefiltered_model_counts
            .splitCsv()
            .filter { row -> dashboard_metric_checkpoint_genomes.contains(row[0]) }
            .map { row -> tuple('prefiltered', row[0], file(row[1])) }

        checkpoint_postfiltered_model_counts_dashboard_metric_ch = LOAD_LORA_CHECKPOINTS.out.postfiltered_model_counts
            .splitCsv()
            .filter { row -> dashboard_metric_checkpoint_genomes.contains(row[0]) }
            .map { row -> tuple('postfiltered', row[0], file(row[1])) }

        checkpoint_prefiltered_sqanti3_counts_dashboard_metric_ch = LOAD_LORA_CHECKPOINTS.out.prefiltered_sqanti3_counts
            .splitCsv()
            .filter { row -> dashboard_metric_checkpoint_genomes.contains(row[0]) }
            .map { row -> tuple('prefiltered', row[0], file(row[1])) }

        checkpoint_postfiltered_sqanti3_counts_dashboard_metric_ch = LOAD_LORA_CHECKPOINTS.out.postfiltered_sqanti3_counts
            .splitCsv()
            .filter { row -> dashboard_metric_checkpoint_genomes.contains(row[0]) }
            .map { row -> tuple('postfiltered', row[0], file(row[1])) }

        checkpoint_full_length_stats_dashboard_metric_ch = LOAD_LORA_CHECKPOINTS.out.full_length_stats
            .splitCsv()
            .filter { row -> dashboard_metric_checkpoint_genomes.contains(row[0]) }
            .map { row -> tuple(row[0], file(row[1])) }

        checkpoint_data_checks_dashboard_ch = LOAD_LORA_CHECKPOINTS.out.data_checks
            .splitCsv()
            .filter { row -> dashboard_metric_checkpoint_genomes.contains(row[0]) }
            .map { row -> tuple(row[0], file(row[1])) }

        checkpoint_dashboard_bam_stats_inputs_ch = LOAD_LORA_CHECKPOINTS.out.dashboard_bam_stats_inputs
            .splitCsv()
            .filter { row -> !dashboard_metric_checkpoint_genomes.contains(row[1]) }
            .map { row -> tuple(row[0], row[1], row[2], file(row[3]), file(row[4])) }

        read_id_maps_ch = Channel.empty()

        MAPPING_READS_CHECKPOINT(checkpoint_read_prep_fasta_ch)

        read_prep_bam_ch = PREPARING_PIPELINE.out.bam_ch
            .filter { reads, reference, genome_id, raw_file_count -> read_prep_genomes.contains(genome_id) }

        SAMTOOLS_CHECKPOINT(MAPPING_READS_CHECKPOINT.out.mapped_reads, read_prep_bam_ch)

        SAMTOOLS_FROM_MERGED_BAM(checkpoint_mapped_bam_ch)

        bam_to_gtf_prepared_bam_ch = SAMTOOLS_CHECKPOINT.out.prepared_bam.mix(SAMTOOLS_FROM_MERGED_BAM.out.prepared_bam)
        bam_to_gtf_chromosome_bam_ch = SAMTOOLS_CHECKPOINT.out.chromosome_bam.mix(SAMTOOLS_FROM_MERGED_BAM.out.chromosome_bam)

        BAM_TO_GTF_CHECKPOINT(bam_to_gtf_prepared_bam_ch, bam_to_gtf_chromosome_bam_ch)

        MERGING_TRANSCRIPTS_CHECKPOINT(BAM_TO_GTF_CHECKPOINT.out.prepared_gtf)

        full_length_chromosome_bam_ch = BAM_TO_GTF_CHECKPOINT.out.chromosome_bam.mix(checkpoint_canonical_bam_ch)
        full_length_chromosome_models_ch = MERGING_TRANSCRIPTS_CHECKPOINT.out.transcript_models_per_chromosome.mix(checkpoint_chromosome_models_ch)

        FULL_LENGTH_FILTER_CHECKPOINT(full_length_chromosome_bam_ch, full_length_chromosome_models_ch)

        prefiltered_models_ch = MERGING_TRANSCRIPTS_CHECKPOINT.out.prefiltered_models.mix(checkpoint_prefiltered_models_ch)
        prefiltered_contains_ch = MERGING_TRANSCRIPTS_CHECKPOINT.out.prefiltered_contains.mix(checkpoint_prefiltered_contains_ch)
        prefiltered_models_per_chromosome_ch = MERGING_TRANSCRIPTS_CHECKPOINT.out.transcript_models_per_chromosome.mix(checkpoint_chromosome_models_ch)
        filtered_models_ch = FULL_LENGTH_FILTER_CHECKPOINT.out.filtered_models.mix(checkpoint_filtered_models_ch)
        filtered_models_per_chromosome_ch = FULL_LENGTH_FILTER_CHECKPOINT.out.filtered_models_per_chromosome.mix(checkpoint_filtered_chromosome_models_ch)
        full_length_stats_ch = FULL_LENGTH_FILTER_CHECKPOINT.out.full_length_stats.mix(checkpoint_full_length_stats_ch)

        dashboard_prefiltered_models_ch = MERGING_TRANSCRIPTS_CHECKPOINT.out.prefiltered_models.mix(checkpoint_prefiltered_models_dashboard_ch)
        dashboard_prefiltered_contains_ch = MERGING_TRANSCRIPTS_CHECKPOINT.out.prefiltered_contains.mix(checkpoint_prefiltered_contains_dashboard_ch)
        dashboard_filtered_models_ch = FULL_LENGTH_FILTER_CHECKPOINT.out.filtered_models.mix(checkpoint_filtered_models_dashboard_ch)
        dashboard_full_length_stats_ch = FULL_LENGTH_FILTER_CHECKPOINT.out.full_length_stats.mix(checkpoint_full_length_stats_dashboard_ch)
        dashboard_read_summaries_ch = checkpoint_read_summaries_dashboard_ch
        rerun_bam_stats_inputs_ch = SAMTOOLS_CHECKPOINT.out.sorted_bam_with_index
            .map { genome_id, bam, bai, _reference, _raw_file_count -> tuple('sample', genome_id, bam.simpleName.toString(), bam, bai) }
            .mix(SAMTOOLS_CHECKPOINT.out.merged_bam
                .map { bam, bai, _reference, genome_id -> tuple('merged', genome_id, genome_id, bam, bai) })
            .mix(SAMTOOLS_CHECKPOINT.out.chromosome_bam
                .map { bam, bai, _reference, genome_id, chromosome_id, _chromosome_count, _bin_count -> tuple('chromosome', genome_id, chromosome_id, bam, bai) })
            .mix(SAMTOOLS_FROM_MERGED_BAM.out.chromosome_bam
                .map { bam, bai, _reference, genome_id, chromosome_id, _chromosome_count, _bin_count -> tuple('chromosome', genome_id, chromosome_id, bam, bai) })
            .mix(BAM_TO_GTF_CHECKPOINT.out.chromosome_bam
                .map { bam, bai, genome_id, chromosome_id, _chromosome_count -> tuple('canonical', genome_id, chromosome_id, bam, bai) })
        dashboard_bam_stats_inputs_ch = checkpoint_dashboard_bam_stats_inputs_ch.mix(rerun_bam_stats_inputs_ch)
    }

    if (checkpoint_resume == 'none') {
        sqanti_prefiltered_models_ch = prefiltered_models_ch
        sqanti_filtered_models_ch = filtered_models_ch
        sqanti_prefiltered_models_per_chromosome_ch = prefiltered_models_per_chromosome_ch
        sqanti_filtered_models_per_chromosome_ch = filtered_models_per_chromosome_ch
    } else {
        sqanti_prefiltered_models_ch = MERGING_TRANSCRIPTS_CHECKPOINT.out.prefiltered_models
            .filter { genome_id, _prefiltered_gtf -> checkpoint_plan[genome_id] in ['read_prep', 'mapped_bam'] }
            .mix(checkpoint_prefiltered_models_ch
                .filter { genome_id, _prefiltered_gtf -> checkpoint_plan[genome_id] in ['merged_models', 'full_length'] })

        sqanti_filtered_models_ch = FULL_LENGTH_FILTER_CHECKPOINT.out.filtered_models
            .filter { genome_id, _postfiltered_gtf, _postfiltered_contains -> checkpoint_plan[genome_id] in ['read_prep', 'mapped_bam', 'merged_models'] }
            .mix(checkpoint_filtered_models_ch
                .filter { genome_id, _postfiltered_gtf, _postfiltered_contains -> checkpoint_plan[genome_id] == 'full_length' })

        sqanti_prefiltered_models_per_chromosome_ch = MERGING_TRANSCRIPTS_CHECKPOINT.out.transcript_models_per_chromosome
            .filter { genome_id, _gtf, _contains, _chromosome_count, _chromosome_id -> checkpoint_plan[genome_id] in ['read_prep', 'mapped_bam'] }
            .mix(checkpoint_chromosome_models_ch
                .filter { genome_id, _gtf, _contains, _chromosome_count, _chromosome_id -> checkpoint_plan[genome_id] in ['merged_models', 'full_length'] })

        sqanti_filtered_models_per_chromosome_ch = FULL_LENGTH_FILTER_CHECKPOINT.out.filtered_models_per_chromosome
            .filter { genome_id, _chromosome_count, _chromosome_id, _gtf, _contains -> checkpoint_plan[genome_id] in ['read_prep', 'mapped_bam', 'merged_models'] }
            .mix(checkpoint_filtered_chromosome_models_ch
                .filter { genome_id, _chromosome_count, _chromosome_id, _gtf, _contains -> checkpoint_plan[genome_id] == 'full_length' })
    }

    run_sqanti3_workflow = run_sqanti3 && (checkpoint_resume == 'none' || sqanti_run_genomes)
    if (run_sqanti3_workflow) {
        SQANTI3(
            PREPARING_PIPELINE.out.sqanti3_reference_ch,
            sqanti_prefiltered_models_ch,
            sqanti_filtered_models_ch,
            sqanti_prefiltered_models_per_chromosome_ch,
            sqanti_filtered_models_per_chromosome_ch
        )
    }

    if (run_dashboard) {
        if (run_sqanti3_workflow) {
            sqanti3_prefiltered_results = SQANTI3.out.sqanti3_results
                .filter { stage, genome_id, sqanti3_dir -> stage == 'prefiltered' }
                .map { stage, genome_id, sqanti3_dir -> tuple(genome_id, sqanti3_dir) }

            sqanti3_postfiltered_results = SQANTI3.out.sqanti3_results
                .filter { stage, genome_id, sqanti3_dir -> stage == 'postfiltered' }
                .map { stage, genome_id, sqanti3_dir -> tuple(genome_id, sqanti3_dir) }

            sqanti3_dashboard_inputs = sqanti3_prefiltered_results
                .join(sqanti3_postfiltered_results, by: 0)
                .map { genome_id, prefiltered_sqanti3, postfiltered_sqanti3 ->
                    tuple(genome_id, prefiltered_sqanti3, postfiltered_sqanti3)
            }
            PREPARE_SQANTI3_DASHBOARD_INPUT(sqanti3_dashboard_inputs)
            sqanti3_dashboard_ch = PREPARE_SQANTI3_DASHBOARD_INPUT.out.sqanti3_dashboard_input.mix(sqanti3_existing_dashboard_ch)
        } else {
            PREPARE_EMPTY_SQANTI3_DASHBOARD_INPUT(
                dashboard_empty_sqanti3_genomes_ch
                    .filter { genome_id -> !sqanti_existing_genomes.contains(genome_id) }
            )
            empty_sqanti3_dashboard_ch = PREPARE_EMPTY_SQANTI3_DASHBOARD_INPUT.out.sqanti3_dashboard_input
            sqanti3_dashboard_ch = empty_sqanti3_dashboard_ch.mix(sqanti3_existing_dashboard_ch)
        }

        sample_count_inputs = dashboard_genomes_ch
            .filter { genome_id -> !dashboard_metric_checkpoint_genomes.contains(genome_id) }
            .combine(PREPARING_PIPELINE.out.sample_metadata)
            .map { genome_id, sample_metadata -> tuple(genome_id, sample_metadata) }
        COUNT_LORA_SAMPLES(sample_count_inputs)

        input_read_count_inputs = dashboard_read_summaries_ch
            .filter { genome_id, _read_summary -> !dashboard_metric_checkpoint_genomes.contains(genome_id) }
            .map { genome_id, read_summary ->
                tuple(genome_id, read_summary)
            }
        COUNT_LORA_INPUT_READS(input_read_count_inputs)

        prefiltered_model_count_inputs = dashboard_prefiltered_models_ch
            .filter { genome_id, _prefiltered_gtf -> !dashboard_metric_checkpoint_genomes.contains(genome_id) }
            .join(dashboard_prefiltered_contains_ch, by: 0)
            .map { genome_id, prefiltered_gtf, prefiltered_contains ->
                tuple('prefiltered', genome_id, prefiltered_gtf, prefiltered_contains)
            }

        postfiltered_model_count_inputs = dashboard_filtered_models_ch
            .filter { genome_id, _postfiltered_gtf, _postfiltered_contains -> !dashboard_metric_checkpoint_genomes.contains(genome_id) }
            .map { genome_id, postfiltered_gtf, postfiltered_contains ->
                tuple('postfiltered', genome_id, postfiltered_gtf, postfiltered_contains)
            }

        COUNT_LORA_MODELS(prefiltered_model_count_inputs.mix(postfiltered_model_count_inputs))

        sqanti3_count_inputs = sqanti3_dashboard_ch
            .filter { genome_id, _sqanti3_available, _sqanti3_dir -> !dashboard_metric_checkpoint_genomes.contains(genome_id) }
            .flatMap { genome_id, sqanti3_available, sqanti3_dir ->
                [
                    tuple('prefiltered', genome_id, sqanti3_available, sqanti3_dir),
                    tuple('postfiltered', genome_id, sqanti3_available, sqanti3_dir)
                ]
            }

        COUNT_LORA_SQANTI3(sqanti3_count_inputs)

        sample_counts_for_dashboard = COUNT_LORA_SAMPLES.out.sample_counts
            .map { genome_id, sample_counts -> sample_counts }
            .mix(checkpoint_sample_counts_dashboard_ch
                .map { genome_id, sample_counts -> sample_counts })

        input_read_counts_for_dashboard = COUNT_LORA_INPUT_READS.out.input_read_counts
            .map { genome_id, input_read_counts -> input_read_counts }
            .mix(checkpoint_input_read_counts_dashboard_ch
                .map { genome_id, input_read_counts -> input_read_counts })

        model_counts_for_dashboard = COUNT_LORA_MODELS.out.model_counts
            .map { stage, genome_id, model_counts -> model_counts }
            .mix(checkpoint_prefiltered_model_counts_dashboard_metric_ch
                .map { stage, genome_id, model_counts -> model_counts })
            .mix(checkpoint_postfiltered_model_counts_dashboard_metric_ch
                .map { stage, genome_id, model_counts -> model_counts })

        sqanti3_counts_for_dashboard = COUNT_LORA_SQANTI3.out.sqanti3_counts
            .map { stage, genome_id, sqanti3_counts -> sqanti3_counts }
            .mix(checkpoint_prefiltered_sqanti3_counts_dashboard_metric_ch
                .map { stage, genome_id, sqanti3_counts -> sqanti3_counts })
            .mix(checkpoint_postfiltered_sqanti3_counts_dashboard_metric_ch
                .map { stage, genome_id, sqanti3_counts -> sqanti3_counts })

        full_length_stats_for_dashboard = dashboard_full_length_stats_ch
            .filter { genome_id, _full_length_stats -> !dashboard_metric_checkpoint_genomes.contains(genome_id) }
            .map { genome_id, full_length_stats -> full_length_stats }
            .mix(checkpoint_full_length_stats_dashboard_metric_ch
                .map { genome_id, full_length_stats -> full_length_stats })

        sample_counts_files_for_dashboard = sample_counts_for_dashboard.collect().ifEmpty([])
        input_read_counts_files_for_dashboard = input_read_counts_for_dashboard.collect().ifEmpty([])
        model_counts_files_for_dashboard = model_counts_for_dashboard.collect().ifEmpty([])
        sqanti3_counts_files_for_dashboard = sqanti3_counts_for_dashboard.collect().ifEmpty([])
        full_length_stats_files_for_dashboard = full_length_stats_for_dashboard.collect().ifEmpty([])

        if (dashboard_only_resume) {
            data_checks_for_dashboard = checkpoint_data_checks_dashboard_ch
                .map { genome_id, data_checks -> data_checks }
        } else {
            CHECK_LORA_BAM_STATS(dashboard_bam_stats_inputs_ch)
            bam_stats_files_for_checks = CHECK_LORA_BAM_STATS.out.bam_stats
                .map { genome_id, bam_stats -> bam_stats }
                .collect()
                .ifEmpty([])

            BUILD_LORA_DATA_CHECKS(
                PREPARING_PIPELINE.out.sample_metadata,
                sample_counts_files_for_dashboard,
                input_read_counts_files_for_dashboard,
                model_counts_files_for_dashboard,
                sqanti3_counts_files_for_dashboard,
                full_length_stats_files_for_dashboard,
                bam_stats_files_for_checks
            )

            data_checks_for_dashboard = BUILD_LORA_DATA_CHECKS.out.data_checks
                .flatten()
        }

        BUILD_LORA_DASHBOARD_MANIFESTS(
            sample_counts_files_for_dashboard,
            input_read_counts_files_for_dashboard,
            model_counts_files_for_dashboard,
            sqanti3_counts_files_for_dashboard,
            full_length_stats_files_for_dashboard,
            data_checks_for_dashboard.collect().ifEmpty([]),
            PREPARING_PIPELINE.out.sample_metadata,
            file(reference_sheet)
        )

        dashboard_input_manifest_inputs = BUILD_LORA_DASHBOARD_MANIFESTS.out.dashboard_inputs
            .flatten()
            .map { dashboard_input_dir ->
                tuple(dashboard_input_dir.name.toString(), dashboard_input_dir)
            }

        SUMMARISE_LORA_DASHBOARD(dashboard_input_manifest_inputs)

        dashboard_render_inputs = SUMMARISE_LORA_DASHBOARD.out.dashboard_summary
            .map { genome_id, _run_summary, dashboard_data ->
                tuple(genome_id, dashboard_data)
            }

        BUILD_LORA_DASHBOARD(dashboard_render_inputs)
        BUILD_LORA_DASHBOARD_INDEX(
            BUILD_LORA_DASHBOARD.out.dashboard
                .map { genome_id, _index_html, _dashboard_data -> genome_id }
                .collect()
        )
    }

    emit:
    prefiltered_models = prefiltered_models_ch
    filtered_models = filtered_models_ch
    full_length_stats = full_length_stats_ch
    read_id_maps = read_id_maps_ch
    read_summaries = read_summaries_ch
    sample_metadata = PREPARING_PIPELINE.out.sample_metadata
}

workflow LORA_MAIN {

    take:
    sample_sheet
    reference_sheet

    main:
    LORA(sample_sheet, reference_sheet)

    emit:
    prefiltered_models = LORA.out.prefiltered_models
    filtered_models = LORA.out.filtered_models
    full_length_stats = LORA.out.full_length_stats
    read_id_maps = LORA.out.read_id_maps
    read_summaries = LORA.out.read_summaries
    sample_metadata = LORA.out.sample_metadata
}
