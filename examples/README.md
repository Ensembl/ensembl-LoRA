# LoRA Example Input Sheets

These files are templates. Edit the paths before running LoRA.

LoRA input sheets may include a header row. Headered sheets are recommended because they make optional reporting metadata easier to check. Older headerless sheets still work.

## Sample Sheet

`sample_sheet.csv` uses these columns:

```text
reads,format,seq_source,map_key,species
```

- `reads`: path to a FASTQ, FASTA, or BAM file.
- `format`: `fastq`, `fasta`, or `bam`.
- `seq_source`: `nanopore` or `pacbio`.
- `map_key`: key used to match reads to a reference genome row.
- `species`: species label used to match reads to a reference genome row.

`map_key` and `species` must match a row in the reference sheet.

`sample_sheet_with_metadata.csv` adds optional reporting metadata:

```text
reads,format,seq_source,map_key,species,tissue,cell_line,sample_id
```

`tissue` and `cell_line` may be blank or set to `NA` when they do not apply. `sample_id` controls dashboard sample counts. If it is blank or absent, LoRA infers a sample group from the filename and collapses split replicate parts such as `rep1a`, `rep1b`, and `rep1aaz` to `rep1`.

`sample_sheet_with_read_id_map.csv` adds an optional stable-renaming map:

```text
reads,format,seq_source,map_key,species,tissue,cell_line,sample_id,read_id_map
```

If `read_id_map` is blank or `NA`, LoRA renames FASTA/FASTQ reads sequentially and writes a new compressed compact two-column map: `original_read_id` and `new_read_id`. If `read_id_map` points to a previous per-input `*.read_id_map.tsv.gz` or `*.read_id_map.tsv`, LoRA reuses those original-to-new read ID pairings and writes the effective map for the current run. Older five-column LoRA maps are still accepted as replay inputs. Use `NA` for BAM rows because LoRA does not rename BAM reads in this step.

LoRA also writes per-input `*.read_summary.tsv` files during renaming. These hold input read counts plus method/tissue/cell-line metadata and are used by the dashboard. Combined per-genome read ID maps are disabled by default; set `--combine_read_id_maps true` if you need one.

Headerless sample sheets are still accepted with 5, 7, 8, 9, or 10 columns. Headered sheets are the clearest option when using optional metadata.

## Reference Sheet Without SQANTI3

`reference_sheet_no_sqanti3.csv` uses these columns:

```text
reference_file,genome_id,map_key,species
```

When no `reference_gtf` value is provided, LoRA runs normally and skips SQANTI3 for that genome.

## Reference Sheet With SQANTI3

`reference_sheet_sqanti3.csv` uses these columns:

```text
reference_file,genome_id,map_key,species,reference_gtf,polya_peak,cage_peak
```

The `reference_gtf` column enables SQANTI3 for that reference genome. It must point to a GTF file, not GFF3. Compressed `.gz` FASTA and GTF files are accepted by LoRA's SQANTI3 reference-preparation step.

`polya_peak` and `cage_peak` are optional BED files passed to SQANTI3 as `--polyA_peak` and `--CAGE_peak`. They may be supplied together, individually, or left blank/`NA`, but they require `reference_gtf`.

Valid SQANTI3 evidence combinations:

```text
reference_gtf + polya_peak + cage_peak
reference_gtf + polya_peak
reference_gtf + cage_peak
reference_gtf only
no reference_gtf and no polya_peak/cage_peak
```

SQANTI3 is skipped if `reference_gtf` is absent or set to an empty/placeholder value such as `NA`. LoRA fails early if `polya_peak` or `cage_peak` is provided without `reference_gtf`.

polyA and CAGE BED files must have chromosome in column 1, zero-based start in column 2, one-based end in column 3, and strand in column 6.

Headerless reference sheets are still accepted with 4, 5, 6, or 7 columns. For headerless CAGE-only input, use `NA` in the `polya_peak` position.
