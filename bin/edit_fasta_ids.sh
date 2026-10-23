#!/bin/bash
set -euo pipefail

fasta="$1"
read_id_map="${2:-NA}"
source_file="${3:-$fasta}"
seq_source="${4:-NA}"
genome_id="${5:-NA}"
tissue="${6:-NA}"
cell_line="${7:-NA}"
cpus="${8:-1}"
compression_level="${9:-1}"
fasta_name=$(basename "$fasta")
accession_id=${fasta_name/_fasta.gz/}

normalise_accession_id() {
    local id="$1"
    echo "$id" | sed -E 's/([0-9]+)reads$/\1/'
}

accession_id=$(normalise_accession_id "$accession_id")
map_file="${accession_id}.fasta.read_id_map.tsv"
summary_file="${accession_id}.fasta.read_summary.tsv"
renamed_fasta="${accession_id}.reads.fasta"

compress_file() {
    local file="$1"
    if ! command -v pigz >/dev/null 2>&1; then
        echo "ERROR: pigz is required for LoRA read-prep compression but was not found in PATH" >&2
        exit 1
    fi
    pigz -f -p "$cpus" "-${compression_level}" "$file"
}

rename_reads_with_optional_map.py \
    --input "$fasta" \
    --format fasta \
    --accession_id "$accession_id" \
    --output "$renamed_fasta" \
    --output_map "$map_file" \
    --source_file "$source_file" \
    --summary "$summary_file" \
    --seq_source "$seq_source" \
    --genome_id "$genome_id" \
    --tissue "$tissue" \
    --cell_line "$cell_line" \
    --provided_map "$read_id_map"

compress_file "$renamed_fasta"
compress_file "$map_file"
