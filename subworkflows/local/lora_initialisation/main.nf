//
// LoRA pipeline initialisation: nf-schema help, parameter summary, and validation.
//

include { paramsSummaryLog   } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'
include { paramsHelp         } from 'plugin/nf-schema'

workflow LORA_INITIALISATION {

    take:
    input_workflow
    validate_params
    monochrome_logs
    help
    help_full
    show_hidden
    nextflow_cli_args

    main:
    def before_text = ""
    def after_text = ""
    def command = "nextflow run LoRA -profile slurm,singularity --input samplesheet.csv --reference_sheet references.csv --outdir results --runID run_name"

    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
        after_text = after_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    if (help || help_full) {
        log.info paramsHelp(
            [
                parametersSchema: "${projectDir}/nextflow_schema.json",
                beforeText: before_text,
                afterText: after_text,
                command: command,
                showHidden: show_hidden,
                fullHelp: help_full
            ],
            (params.help instanceof String && params.help != "true") ? params.help : ""
        )
        exit 0
    }

    if (params.version) {
        log.info "${input_workflow.manifest.name} ${input_workflow.manifest.version}"
        exit 0
    }

    if (nextflow_cli_args.any { it == "--outDir" || it.startsWith("--outDir=") }) {
        error "Parameter --outDir is no longer supported. Use --outdir instead."
    }

    if (validate_params) {
        validateParameters([parametersSchema: "${projectDir}/nextflow_schema.json"])
    }

    def resolved_input = params.input ?: params.sample_sheet
    def missing_required = []
    if (!resolved_input) {
        missing_required << "input or sample_sheet"
    }
    if (!params.reference_sheet) {
        missing_required << "reference_sheet"
    }
    if (!params.outdir) {
        missing_required << "outdir"
    }
    if (!params.runID) {
        missing_required << "runID"
    }
    if (missing_required) {
        error "Missing required LoRA parameter(s): ${missing_required.join(', ')}"
    }

    log.info before_text
    log.info paramsSummaryLog([parametersSchema: "${projectDir}/nextflow_schema.json"], input_workflow)
    log.info after_text

    emit:
    dummy_emit = true
}
