# LoRA

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-23aa8c.svg)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.7-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![Run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

A Nextflow pipeline for Long-read RNA-seq Analysis. LoRA processes long-read sequencing data (PacBio and Oxford Nanopore) to generate high-quality transcript annotations through read mapping, alignment processing, transcript assembly, and full-length read filtering.

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Parameters](#parameters)
- [Pipeline Workflow](#pipeline-workflow)
- [Directory Structure](#directory-structure)
- [Output](#output)
- [Development](#development)
- [Configuration](#configuration)
- [Author](#author)
- [License](#license)

## Features

- **Multi-Platform Support**: Processes both PacBio and Oxford Nanopore long-read data
- **Automated Read Processing**: Handles FASTQ to FASTA conversion and read ID standardization
- **Optimized Mapping**: Uses minimap2 with platform-specific settings for accurate alignment
- **Transcript Assembly**: Employs tmerge algorithm for robust transcript model building
- **Full-Length Filtering**: Identifies and filters full-length transcripts based on intron chain matching
- **Chromosome-level Processing**: Per-chromosome analysis for efficient large genome handling
- **Quality Control**: Comprehensive read statistics and filtering metrics
- **Flexible Input**: Supports both raw sequencing data and pre-aligned BAM files

## Requirements

- [Nextflow](https://www.nextflow.io/) (≥23.04.0)
- [Singularity](https://sylabs.io/singularity/) or [Conda](https://docs.conda.io/en/latest/)
- Linux environment (tested on EBI Slurm HPC)

## Quick Start

1. **Prepare input sheets**

   Create a sample sheet (CSV format). Headered sheets are recommended:
   ```csv
   reads,format,seq_source,map_key,species
   /path/to/sample1.fastq.gz,fastq,nanopore,hg38,human
   /path/to/sample2.fasta.gz,fasta,pacbio,hg38,human
   /path/to/sample3.bam,bam,nanopore,hg38,human
   ```

   Sample sheet columns:
   ```text
   reads,format,seq_source,map_key,species
   ```

   You can also add optional reporting metadata:
   ```text
   reads,format,seq_source,map_key,species,tissue,cell_line,sample_id
   ```

   `sample_id` is optional. If it is blank or absent, LoRA infers it from the input filename and collapses split replicate parts such as `rep1a`, `rep1b`, and `rep1aaz` to `rep1` for dashboard sample counting.

   To reproduce previous FASTA/FASTQ read renaming, add an optional per-input read ID map:
   ```text
   reads,format,seq_source,map_key,species,tissue,cell_line,sample_id,read_id_map
   ```

   Leave `read_id_map` blank or `NA` to generate new names normally. Provide a previous per-input `*.read_id_map.tsv.gz` or `*.read_id_map.tsv` to replay the same original-to-LoRA read ID pairings. New maps are compressed compact two-column files with `original_read_id` and `new_read_id`; older five-column LoRA maps are still accepted for replay. BAM rows should use `NA`.

   Create a reference sheet without SQANTI3:
   ```csv
   reference_file,genome_id,map_key,species
   /path/to/genome.fa.gz,hg38,hg38,human
   ```

   Reference sheet columns without SQANTI3:
   ```text
   reference_file,genome_id,map_key,species
   ```

   To enable SQANTI3 for a reference genome, add a `reference_gtf` column:
   ```csv
   reference_file,genome_id,map_key,species,reference_gtf,polya_peak,cage_peak
   /path/to/genome.fa.gz,hg38,hg38,human,/path/to/reference.annotation.gtf.gz,/path/to/polyA_peaks.bed.gz,/path/to/CAGE_peaks.bed.gz
   ```

   Reference sheet columns with SQANTI3:
   ```text
   reference_file,genome_id,map_key,species,reference_gtf,polya_peak,cage_peak
   ```

   `polya_peak` and `cage_peak` are optional BED files for SQANTI3 transcript-end evidence. They can be provided together, individually, or set to `NA`, but LoRA will fail if either is provided without `reference_gtf`.
   BED files should have chromosome in column 1, zero-based start in column 2, one-based end in column 3, and strand in column 6.

   `map_key` and `species` must match between the sample and reference sheets. Older headerless sheets still work. See `examples/` for editable templates.

2. **Configure parameters**
   
   Edit `params.config`:
   ```groovy
   params {
       input = "/path/to/sample_sheet.csv"
       reference_sheet = "/path/to/reference_sheet.csv"
       outdir = "/path/to/output"
       runID = "my_run"
       minFL = 1  // Minimum full-length reads threshold
       genomeSize = "normal"
       species = "human"
   }
   ```

3. **Run the pipeline**
   ```bash
   nextflow run main.nf -profile slurm,singularity
   ```

## Usage

### Basic Usage

```bash
nextflow run main.nf
```

### On Slurm With Singularity

```bash
nextflow run main.nf -profile slurm,singularity
```

### Resume a Previous Run

```bash
nextflow run main.nf -resume
```

### Audit Published Checkpoints

Nextflow `-resume` uses the work directory cache. LoRA can also audit trusted published outputs so you can see how far each genome got:

```bash
bin/lora_checkpoint_status.py --outdir /path/to/output --run_id my_run
```

To restart from a published checkpoint:

```bash
nextflow run main.nf ... --resume_from mapped_bam --resume_scope common
```

Or let LoRA choose the latest common published checkpoint shared by all genomes:

```bash
nextflow run main.nf ... --resume_from auto --resume_scope common
```

### Help

```bash
nextflow run main.nf --help
nextflow run main.nf --help_full
```

## Parameters

| Parameter | Description | Required | Default |
|-----------|-------------|----------|---------|
| `input` / `sample_sheet` | CSV file with read inputs; optional tissue, cell_line, sample_id, and read_id_map columns are supported | Yes | - |
| `reference_sheet` | CSV file with reference inputs; optional reference_gtf column enables SQANTI3, with optional polya_peak/cage_peak BED evidence | Yes | - |
| `sqanti3_by_chromosome` | Run SQANTI3 per chromosome using chromosome-specific model GTFs and reference GTF/BED subsets. Full-stage mode runs one SQANTI3 job per stage. | No | `false` |
| `outdir` | Output directory for all results | Yes | - |
| `runID` | Unique identifier for the pipeline run | Yes | - |
| `minFL` | Minimum number of full-length reads supporting a transcript | No | 1 |
| `genomeSize` | Minimap2 genome mode; `large` enables `--split-prefix=tmp` | No | normal |
| `species` | Species identifier matching reference sheet | Yes | - |
| `combine_read_id_maps` | Combine per-input compact read ID maps into one per-genome map; off by default to avoid very large files | No | false |
| `read_prep_compression_level` | gzip/pigz compression level for renamed read files and compact read ID maps | No | 1 |
| `resume_from` | Published-output checkpoint to resume from: `none`, `auto`, `read_prep`, `mapped_bam`, `bam_to_gtf`, `merged_models`, `full_length`, `sqanti3`, or `dashboard`; `bam_to_gtf` is audit-only and fails with guidance if requested | No | none |
| `resume_scope` | Checkpoint scope; use `common` for the latest shared checkpoint or `per_genome` for mixed per-genome recovery | No | common |
| `validate_params` | Validate parameters with nf-schema | No | true |

See `docs/usage.md` and `docs/output.md` for the modernised Nextflow layout and output guide.

## Pipeline Workflow

The pipeline consists of several sequential stages:

### 1. Pipeline Preparation
- Parses sample and reference sheets
- Validates input files
- Creates processing channels for FASTQ, FASTA, and BAM files

### 2. Read File Preparation
For FASTQ inputs:
- Converts FASTQ directly to renamed FASTA in one pass
- Standardizes read IDs for downstream processing
- Writes compact read ID maps and per-input read summaries
- Calculates read count and average read length during renaming

For FASTA inputs:
- Standardizes read IDs in one pass
- Writes compact read ID maps and per-input read summaries
- Calculates read count and average read length during renaming

### 3. Read Mapping
- **Nanopore reads**: Maps using minimap2 with `-ax splice` preset
- **PacBio reads**: Maps using minimap2 with `-ax splice:hq` preset
- Streams minimap2 output directly into filtered, sorted, indexed BAM to avoid large SAM/intermediate-BAM files
- The minimap2 container must include both `minimap2` and `samtools`; LoRA currently uses `docker://community.wave.seqera.io/library/minimap2_samtools:3a608f92b72cacce`

### 4. BAM Processing with Samtools
- Sorts alignments by coordinate
- Merges multiple BAM files per sample
- Filters for canonical splice sites
- Splits BAM files by chromosome

### 5. BAM to GTF Conversion
- Converts BAM alignments to BED12 format
- Fixes strand information
- Converts BED12 to GenePred format
- Generates GTF files from GenePred
- Splits GTF files by chromosome

### 6. Transcript Merging
First merge (per-chromosome):
- Runs tmerge algorithm on chromosome-specific GTF files
- Combines transcript models across chromosomes
- Extracts intron chain signatures from models

Second merge (final assembly):
- Performs final tmerge across all samples
- Generates comprehensive transcript models
- Outputs final GTF annotation

### 7. Full-Length Read Filtering
- Extracts intron chains from aligned reads
- Matches read intron chains against transcript models
- Identifies full-length reads for each transcript
- Filters transcripts based on `minFL` threshold
- Generates statistics on full-length read support

## Directory Structure

```
LoRA/
├── main.nf                 # Main pipeline workflow
├── nextflow.config         # Nextflow configuration
├── params.config           # Pipeline parameters
├── modules/                # Individual process modules
│   └── local/
│       ├── edit_*_read_ids/          # Read ID standardization and FASTA preparation
│       ├── *_minimap2/               # Platform-specific mapping
│       ├── samtools_*/               # BAM processing
│       ├── bam_to_bed/               # Format conversion
│       ├── *_tmerge*/                # Transcript merging
│       ├── intron_chain_extractor_*/ # Intron chain analysis
│       └── filter_*/                 # Full-length filtering
├── subworkflows/          # Reusable subworkflows
│   └── local/
│       ├── preparing_pipeline/       # Input preparation
│       ├── prepare_fast_files/       # Read file processing
│       ├── mapping/                  # Read alignment
│       ├── samtools/                 # BAM processing
│       ├── bam_to_gtf/              # GTF generation
│       ├── merging_transcripts/      # Transcript assembly
│       └── full_length_filter/       # FL filtering
├── bin/                   # Helper scripts
│   ├── process_sample_sheets.py     # Input parsing
│   ├── tmerge.pl                    # Transcript merging
│   ├── tmerge_to_gtf.py             # Format conversion
│   ├── extract_intron_chains*.py    # Intron analysis
│   └── full_length_read_extraction.py # FL filtering
├── conf/                  # Configuration files
│   └── ebi_slurm.config   # EBI HPC Slurm configuration
└── LICENSE                # License information
```

## Output

The pipeline generates:

- **Mapped Reads**: Sorted, indexed BAM files with alignment statistics
- **GTF Annotations**: 
  - Per-chromosome GTF files
  - Merged transcript models
  - Full-length filtered annotations
- **Intron Chain Files**: Intron signatures for reads and transcript models
- **Full-Length Statistics**:
  - Per-transcript FL read counts
  - Filtering summary statistics
  - Contains files mapping reads to transcripts
- **Quality Metrics**:
  - Average read lengths
  - Mapping statistics
  - Transcript support levels
  - Reusable dashboard metric TSVs
  - Data-check TSVs grouped by pipeline stage
- **Dashboard**:
  - Overview, Data Checks, Models, SQANTI3, Files, and Help tabs

### Output Directory Structure

```
outdir/
├── <genome_id>/
│   ├── reads/             # Renamed reads, read summaries, read ID maps, BAMs
│   ├── all_models/        # Prefiltered merged transcript models
│   ├── full_length/       # FL-filtered models, contains files, and stats
│   ├── sqanti3/           # Optional SQANTI3 raw outputs
│   ├── summary/           # Reusable dashboard metric and data-check TSVs
│   └── dashboard/         # Per-genome dashboard
└── dashboard/             # Cross-genome dashboard index
```

## Development

### Module Structure

All modules follow the standard Nextflow DSL2 structure:
```
modules/local/<module_name>/
└── main.nf
```

### Subworkflow Structure

Subworkflows combine multiple modules:
```
subworkflows/local/<subworkflow_name>/
└── main.nf
```

## Configuration

### Profiles

- **slurm**: EBI HPC Slurm executor and resource configuration.
- **singularity**: Singularity container runtime.
- **apptainer**: Apptainer container runtime.
- **conda**: Conda-based software environment.
- **test**: Minimal local test configuration.

### Execution Options

- Local execution is the default when no executor profile is supplied.
- HPC execution with Singularity uses `-profile slurm,singularity`.
- Conda environments are enabled only with `-profile conda`.
- Slurm runs use lenient process caching for robust shared-filesystem execution.
- Retry strategy: Automatic retry with increased resources

## Author

Ryan Merritt

## License

See [LICENSE](LICENSE) file for details.

## Acknowledgments

Built with [Nextflow](https://www.nextflow.io/) - A DSL for data-driven computational pipelines.

### Key Tools

- [minimap2](https://github.com/lh3/minimap2) - Long-read alignment
- [samtools](http://www.htslib.org/) - SAM/BAM processing
- [bedtools](https://bedtools.readthedocs.io/) - Genome arithmetic
- tmerge - Transcript merging algorithm
