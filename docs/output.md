# LoRA Output

LoRA publishes outputs under `params.outdir`.

Main output groups:

- `<genome_id>/reads/`: merged and chromosome-split read alignments.
- `<genome_id>/all_models/`: post-`tmerge` prefiltered transcript models and contains files.
- `<genome_id>/full_length/`: full-length-filtered models, reads, contains files, and stats.
- `<genome_id>/sqanti3/`: optional SQANTI3 reference prep, prefiltered QC, and postfiltered QC outputs.
- `<genome_id>/summary/metrics/`: reusable metric TSVs for samples, input reads, models, SQANTI3, and BAM checks.
- `<genome_id>/summary/checks/`: reusable data-check TSVs with `stage`, `check_name`, `status`, `observed`, `expected`, `message`, and `source_file`.
- `<genome_id>/dashboard/`: per-genome static LoRA dashboard.
- `dashboard/`: top-level dashboard index.
- `pipeline_info/`: Nextflow timeline, report, trace, DAG, and version metadata where available.

SQANTI3 is skipped for any reference genome without a `reference_gtf` value in the reference sheet. Optional `polya_peak` and `cage_peak` BED files are passed to SQANTI3 only when supplied with `reference_gtf`.

If `--sqanti3_by_chromosome true`, LoRA runs SQANTI3 per chromosome, publishes merged stage-level SQANTI3 files in the usual `prefiltered/` and `postfiltered/` directories, and keeps the per-chromosome SQANTI3 outputs under `chunks/<chromosome_id>/` inside each stage directory. Per-chromosome SQANTI3 HTML/PDF reports are skipped with `--report skip` to avoid report-rendering failures after successful classification.

The dashboard includes a `Data Checks` tab built only from precomputed check TSVs. These checks are diagnostic: `warning` marks suspicious but non-fatal observations, while `failed` marks outputs or counts that should be investigated before trusting that stage.
