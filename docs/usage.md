# LoRA Usage

LoRA uses CSV input sheets. Headered sheets are recommended, and older headerless sheets are still accepted.

## Sample Sheet

Core columns:

```text
reads,format,seq_source,map_key,species
```

Optional reporting metadata columns:

```text
tissue,cell_line,sample_id
```

Optional stable read-renaming column:

```text
read_id_map
```

Allowed values:

- `format`: `fastq`, `fasta`, or `bam`
- `seq_source`: `nanopore` or `pacbio`

`map_key` and `species` must match a row in the reference sheet.

`tissue` and `cell_line` are optional. Use blank, `NA`, `none`, `null`, or `.` if they do not apply.

`sample_id` is optional and controls dashboard sample counting. If it is absent or blank, LoRA infers it from the read filename. Split replicate suffixes are collapsed by removing letters immediately after `repN`, so `rep1a`, `rep1b`, `rep1aa`, and `rep1aaz` all count as `rep1`. Physical file counts are still reported separately as input files.

`read_id_map` is optional for FASTA/FASTQ rows. If absent, LoRA keeps its normal sequential renaming and creates a new compressed compact `*.read_id_map.tsv.gz` with `original_read_id` and `new_read_id`. If supplied, it should point to a previous per-input LoRA read ID map, either `.tsv.gz` or `.tsv`, and LoRA will reuse those original-to-new pairings. Current compact maps are matched by `original_read_id`; older five-column maps are matched by `original_header` first and then `original_read_id`. Use `NA` for BAM rows.

LoRA also writes per-input `*.read_summary.tsv` files while renaming reads. These summaries contain `source_file`, `format`, `seq_source`, `genome_id`, `tissue`, `cell_line`, `read_count`, and `average_read_length`, and the dashboard uses them instead of scanning large read ID maps.

Set `--combine_read_id_maps true` only when you need one combined per-genome map. The default is `false` because large runs can still create very large combined maps even after compacting the columns. Per-input compact maps are still produced and can be supplied back through `read_id_map` for stable read-name replay.

Read preparation requires `pigz` and uses compression level `1` by default for speed. Set `--read_prep_compression_level` from `1` to `9` if you want smaller renamed-read/map outputs at the cost of longer read-prep jobs.

Headerless sample sheets are still accepted with 5, 7, 8, 9, or 10 columns, but headered sheets are recommended for optional metadata.

## Reference Sheet

Without SQANTI3:

```text
reference_file,genome_id,map_key,species
```

With SQANTI3:

```text
reference_file,genome_id,map_key,species,reference_gtf,polya_peak,cage_peak
```

The `reference_gtf` column must be a GTF file. Compressed `.gz` FASTA, GTF, polyA BED, and CAGE BED files are accepted by the SQANTI3 reference-preparation step.

`polya_peak` and `cage_peak` are optional BED files passed to SQANTI3 as `--polyA_peak` and `--CAGE_peak`. Either or both can be supplied, but LoRA will fail early if they are provided without `reference_gtf`.

For large model sets, `--sqanti3_by_chromosome true` runs SQANTI3 per chromosome using LoRA's chromosome-level model GTFs. LoRA also filters the SQANTI3 reference GTF and optional polyA/CAGE BED files to the same chromosome before each SQANTI3 task. The merged classification and junction files keep the usual stage-level filenames for downstream LoRA reports. Chromosome SQANTI3 uses `--report skip`; native SQANTI3 HTML/PDF reports are only generated when `--sqanti3_by_chromosome false`.

The BED files should have chromosome in column 1, zero-based start in column 2, one-based end in column 3, and strand in column 6.

Headerless reference sheets are still accepted with 4, 5, 6, or 7 columns. For headerless CAGE-only input, use `NA` in the `polya_peak` position.

## Mapping Container

LoRA streams minimap2 output directly through `samtools view` and `samtools sort` during mapping. The minimap2 container must therefore include both `minimap2` and `samtools`; LoRA currently uses `docker://community.wave.seqera.io/library/minimap2_samtools:3a608f92b72cacce`. Missing tools should be fixed in the container rather than handled by a fallback pipeline path.

## Run

Run LoRA:

```bash
nextflow run main.nf \
  -profile slurm,singularity \
  --input examples/sample_sheet.csv \
  --reference_sheet examples/reference_sheet_sqanti3.csv \
  --outdir results \
  --runID test \
  --genomeSize normal
```

Show schema-driven help:

```bash
nextflow run main.nf --help
nextflow run main.nf --help_full
```

## Published Checkpoint Audit And Resume

Nextflow `-resume` is still available, but it depends on the work directory cache. To see how far a published LoRA output directory got, run:

```bash
bin/lora_checkpoint_status.py --outdir results --run_id test
```

The audit reports each genome across:

```text
read_prep, mapped_bam, bam_to_gtf, merged_models, full_length, sqanti3, dashboard
```

Each stage is labelled `complete`, `partial`, `missing`, `skipped`, or `not_resumable`, and the report includes both a common recommendation across all genomes and a per-genome recommendation. `bam_to_gtf` is audit-only because LoRA does not publish the exact fixed per-bin GTF channel needed to restart there.

Resume from the latest checkpoint shared by all genomes:

```bash
nextflow run main.nf \
  -profile slurm,singularity \
  --input examples/sample_sheet.csv \
  --reference_sheet examples/reference_sheet_sqanti3.csv \
  --outdir results \
  --runID test \
  --resume_from auto \
  --resume_scope common
```

For mixed multi-genome recovery, use:

```bash
nextflow run main.nf ... --resume_from auto --resume_scope per_genome
```

Supported workflow resume checkpoints are `read_prep`, `mapped_bam`, `merged_models`, `full_length`, `sqanti3`, and `dashboard`. Forced `--resume_from bam_to_gtf` fails early with guidance to resume from `mapped_bam` or `merged_models`.
