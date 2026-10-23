#!/usr/bin/env python3
import argparse
import csv
import os
import re
import shutil
from collections import defaultdict


HEADER_MERGE_SUFFIXES = {
    "_classification.txt",
    "_junctions.txt",
    "_isoform_hits.txt",
}

CONCAT_SUFFIXES = {
    "_corrected.gtf",
    "_corrected.fasta",
    "_corrected.cds.gff3",
    "_corrected.faa",
    "_corrected.pep",
}


def count_data_rows(path):
    if not os.path.exists(path):
        return 0
    with open(path) as handle:
        rows = sum(1 for _line in handle)
    return max(0, rows - 1)


def split_sort_key(value):
    name = os.path.basename(str(value))
    match = re.search(r"chunk(\d+)", name)
    if match:
        return (0, int(match.group(1)), name)
    chrom_match = re.fullmatch(r"(?:chr)?(\d+)", name)
    if chrom_match:
        return (1, int(chrom_match.group(1)), name)
    special = {"chrX": 23, "X": 23, "chrY": 24, "Y": 24, "chrM": 25, "MT": 25, "M": 25}
    if name in special:
        return (1, special[name], name)
    parts = [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", name)]
    return (2, parts, name)


def discover_chunk_dirs(root):
    chunk_dirs = []
    for dirpath, _dirnames, filenames in os.walk(root, followlinks=True):
        if any(name.endswith("_classification.txt") for name in filenames):
            chunk_dirs.append(dirpath)
    return sorted(set(chunk_dirs), key=split_sort_key)


def split_id_for_file(filename, base_prefix):
    suffix = suffix_for_file(filename, base_prefix)
    if not suffix:
        return None
    if filename == f"{base_prefix}{suffix}":
        return "full"
    if filename.startswith(f"{base_prefix}_"):
        return filename[len(base_prefix) + 1: -len(suffix)]
    return None


def discover_chunk_files(root, base_prefix):
    chunk_files = defaultdict(list)
    for dirpath, _dirnames, filenames in os.walk(root, followlinks=True):
        for filename in filenames:
            split_id = split_id_for_file(filename, base_prefix)
            if not split_id:
                continue
            chunk_files[split_id].append(os.path.join(dirpath, filename))
    return {chunk_id: sorted(paths, key=split_sort_key) for chunk_id, paths in chunk_files.items()}


def suffix_for_file(filename, base_prefix):
    known_suffixes = sorted(HEADER_MERGE_SUFFIXES | CONCAT_SUFFIXES | {"_merge_manifest.tsv", ".qc_params.txt"}, key=len, reverse=True)
    for suffix in known_suffixes:
        if filename == f"{base_prefix}{suffix}":
            return suffix
        if filename.startswith(f"{base_prefix}_") and filename.endswith(suffix):
            return suffix
    return None


def merge_header_files(paths, output_path):
    wrote_header = False
    with open(output_path, "w") as out:
        for path in sorted(paths, key=split_sort_key):
            with open(path) as handle:
                for line_number, line in enumerate(handle):
                    if line_number == 0:
                        if not wrote_header:
                            out.write(line)
                            wrote_header = True
                    else:
                        out.write(line)


def concatenate_files(paths, output_path):
    with open(output_path, "w") as out:
        for path in sorted(paths, key=split_sort_key):
            with open(path) as handle:
                shutil.copyfileobj(handle, out)


def read_chunk_manifests(manifest_root, expected_genome_id, expected_stage):
    manifests = []
    for dirpath, _dirnames, filenames in os.walk(manifest_root, followlinks=True):
        for filename in filenames:
            if filename.endswith(".sqanti3_chunk_manifest.tsv"):
                manifests.append(os.path.join(dirpath, filename))

    if not manifests:
        raise SystemExit(f"No SQANTI3 chunk manifest files found under {manifest_root}")

    chunks = {}
    for manifest_path in sorted(manifests, key=split_sort_key):
        with open(manifest_path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"stage", "genome_id", "chunk_id", "output_dir"}
            missing = required.difference(reader.fieldnames or [])
            if missing:
                raise SystemExit(
                    f"SQANTI3 chunk manifest {manifest_path} is missing column(s): "
                    f"{', '.join(sorted(missing))}"
                )
            rows = list(reader)
        if len(rows) != 1:
            raise SystemExit(f"SQANTI3 chunk manifest {manifest_path} should contain exactly one data row")

        row = rows[0]
        row_stage = row["stage"]
        row_genome_id = row["genome_id"]
        chunk_id = row["chunk_id"]
        output_dir = row["output_dir"]
        if row_stage != expected_stage or row_genome_id != expected_genome_id:
            raise SystemExit(
                f"SQANTI3 chunk manifest {manifest_path} belongs to "
                f"stage='{row_stage}', genome='{row_genome_id}', but merge expected "
                f"stage='{expected_stage}', genome='{expected_genome_id}'"
            )
        key = (row_stage, row_genome_id, chunk_id)
        if key in chunks:
            raise SystemExit(
                f"Duplicate SQANTI3 chunk manifest for stage='{row_stage}', "
                f"genome='{row_genome_id}', chunk='{chunk_id}': "
                f"{chunks[key]} and {manifest_path}"
            )
        if not os.path.isdir(output_dir):
            raise SystemExit(
                f"SQANTI3 chunk output directory from manifest does not exist or is not a directory: "
                f"{output_dir}"
            )
        chunks[key] = output_dir

    return {key[2]: output_dir for key, output_dir in sorted(chunks.items(), key=lambda item: split_sort_key(item[0][2]))}


def discover_files_in_dirs(chunk_dirs, base_prefix):
    chunk_files = defaultdict(list)
    for chunk_id, chunk_dir in chunk_dirs.items():
        for dirpath, _dirnames, filenames in os.walk(chunk_dir, followlinks=True):
            for filename in filenames:
                if split_id_for_file(filename, base_prefix):
                    chunk_files[chunk_id].append(os.path.join(dirpath, filename))
    return {chunk_id: sorted(paths, key=split_sort_key) for chunk_id, paths in chunk_files.items()}


def merge_outputs(manifest_root, output_dir, run_id, genome_id, stage):
    base_prefix = f"{run_id}_{genome_id}_{stage}"
    os.makedirs(output_dir, exist_ok=True)
    chunk_out = os.path.join(output_dir, "chunks")
    os.makedirs(chunk_out, exist_ok=True)

    chunk_dirs = read_chunk_manifests(manifest_root, genome_id, stage)
    chunk_files = discover_files_in_dirs(chunk_dirs, base_prefix)
    if not chunk_files:
        raise SystemExit(
            f"No SQANTI3 chunk classification files found from manifests under {manifest_root}. "
            "Check the output_dir paths recorded in the *.sqanti3_chunk_manifest.tsv files."
        )

    by_suffix = defaultdict(list)
    manifest_rows = []

    for chunk_name, paths in sorted(chunk_files.items(), key=lambda item: split_sort_key(item[0])):
        copied_dir = os.path.join(chunk_out, chunk_name)
        if os.path.exists(copied_dir):
            shutil.rmtree(copied_dir)
        os.makedirs(copied_dir, exist_ok=True)

        classification_rows = 0
        junction_rows = 0
        source_dirs = set()
        for path in paths:
            if not os.path.isfile(path):
                continue
            filename = os.path.basename(path)
            source_dirs.add(os.path.dirname(path))
            shutil.copyfile(path, os.path.join(copied_dir, filename))
            suffix = suffix_for_file(filename, base_prefix)
            if not suffix:
                continue
            by_suffix[suffix].append(path)
            if suffix == "_classification.txt":
                classification_rows = count_data_rows(path)
            elif suffix == "_junctions.txt":
                junction_rows = count_data_rows(path)
        manifest_rows.append((chunk_name, ",".join(sorted(source_dirs)), classification_rows, junction_rows))

    for suffix, paths in sorted(by_suffix.items()):
        output_path = os.path.join(output_dir, f"{base_prefix}{suffix}")
        if suffix in HEADER_MERGE_SUFFIXES:
            merge_header_files(paths, output_path)
        elif suffix in CONCAT_SUFFIXES:
            concatenate_files(paths, output_path)

    params_files = []
    for chunk_dir in chunk_dirs.values():
        for dirpath, _dirnames, filenames in os.walk(chunk_dir, followlinks=True):
            for filename in filenames:
                if filename.endswith(".qc_params.txt"):
                    params_files.append(os.path.join(dirpath, filename))
    if params_files:
        shutil.copyfile(sorted(params_files, key=split_sort_key)[0], os.path.join(output_dir, f"{base_prefix}.qc_params.txt"))

    manifest_path = os.path.join(output_dir, f"{base_prefix}_merge_manifest.tsv")
    with open(manifest_path, "w") as manifest:
        manifest.write("chunk_id\tchunk_dir\tclassification_rows\tjunction_rows\n")
        for row in sorted(manifest_rows, key=lambda item: item[0]):
            manifest.write("\t".join(str(value) for value in row) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Merge chunked SQANTI3 outputs into stage-level LoRA outputs.")
    parser.add_argument("--manifest_root", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--run_id", required=True)
    parser.add_argument("--genome_id", required=True)
    parser.add_argument("--stage", required=True)
    args = parser.parse_args()
    merge_outputs(args.manifest_root, args.output_dir, args.run_id, args.genome_id, args.stage)


if __name__ == "__main__":
    main()
