#!/usr/bin/env python3
import argparse
import csv
import os


REQUIRED_COLUMNS = [
    "stage",
    "genome_id",
    "chunk_id",
    "isoforms",
    "reference_fasta",
    "reference_gtf",
    "polya_peak",
    "polya_source",
    "cage_peak",
    "cage_source",
]


def main():
    args = get_args()
    rows = read_manifest_rows(args.manifest_root)
    validate_rows(rows)
    write_rows(args.output, rows)


def get_args():
    parser = argparse.ArgumentParser(description="Validate SQANTI3 QC input manifests.")
    parser.add_argument("--manifest_root", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_manifest_rows(manifest_root):
    manifest_paths = []
    for dirpath, _dirnames, filenames in os.walk(manifest_root, followlinks=True):
        for filename in filenames:
            if filename.endswith(".sqanti3_input_manifest.tsv"):
                manifest_paths.append(os.path.join(dirpath, filename))

    if not manifest_paths:
        raise SystemExit(f"No SQANTI3 input manifest files found under {manifest_root}")

    rows = []
    for manifest_path in sorted(manifest_paths):
        with open(manifest_path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = set(REQUIRED_COLUMNS).difference(reader.fieldnames or [])
            if missing:
                raise SystemExit(
                    f"SQANTI3 input manifest {manifest_path} is missing column(s): "
                    f"{', '.join(sorted(missing))}"
                )
            manifest_rows = list(reader)
        if len(manifest_rows) != 1:
            raise SystemExit(f"SQANTI3 input manifest {manifest_path} should contain exactly one data row")
        row = {column: manifest_rows[0].get(column, "") for column in REQUIRED_COLUMNS}
        row["_manifest_path"] = manifest_path
        rows.append(row)
    return rows


def validate_rows(rows):
    task_keys = {}
    for row in rows:
        key = (row["stage"], row["genome_id"], row["chunk_id"])
        if key in task_keys:
            previous = task_keys[key]
            raise SystemExit(
                "Duplicate SQANTI3 input for "
                f"stage='{row['stage']}', genome='{row['genome_id']}', chunk='{row['chunk_id']}'. "
                f"Manifests: {previous['_manifest_path']} and {row['_manifest_path']}. "
                f"Isoforms: {previous['isoforms']} and {row['isoforms']}"
            )
        task_keys[key] = row

        for column in ["isoforms", "reference_fasta", "reference_gtf", "polya_peak", "cage_peak"]:
            value = row[column]
            if value and not os.path.exists(value):
                raise SystemExit(f"SQANTI3 input path does not exist for column '{column}': {value}")


def write_rows(output_path, rows):
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REQUIRED_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in sorted(rows, key=lambda item: (item["stage"], item["genome_id"], item["chunk_id"])):
            writer.writerow({column: row[column] for column in REQUIRED_COLUMNS})


if __name__ == "__main__":
    main()
