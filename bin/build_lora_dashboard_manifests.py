#!/usr/bin/env python3
import argparse
import csv
import os
import shutil
from pathlib import Path


REQUIRED_FIELDS = [
    "sample_counts",
    "input_read_counts",
    "prefiltered_model_counts",
    "postfiltered_model_counts",
    "prefiltered_sqanti3_counts",
    "postfiltered_sqanti3_counts",
    "full_length_stats",
    "data_checks",
]


def main():
    args = get_args()
    output_root = Path(args.output_dir) / "dashboard_inputs"
    output_root.mkdir(parents=True, exist_ok=True)

    by_genome = {}
    add_simple_files(by_genome, args.run_id, args.sample_counts, "sample_counts", "_sample_counts.tsv")
    add_simple_files(by_genome, args.run_id, args.input_read_counts, "input_read_counts", "_input_read_counts.tsv")
    add_simple_files(by_genome, args.run_id, args.full_length_stats, "full_length_stats", "_full_length_stats.tsv")
    add_stage_files(by_genome, args.run_id, args.model_counts, "model_counts.tsv")
    add_stage_files(by_genome, args.run_id, args.sqanti3_counts, "sqanti3_counts.tsv")
    add_simple_files(by_genome, args.run_id, args.data_checks, "data_checks", "_data_checks.tsv")

    if not by_genome:
        raise SystemExit("No dashboard metric files were provided")

    for genome_id in sorted(by_genome):
        genome_files = by_genome[genome_id]
        missing = [field for field in REQUIRED_FIELDS if field not in genome_files]
        if missing:
            raise SystemExit(f"Missing dashboard metric for {genome_id}: {', '.join(missing)}")

        bundle_dir = output_root / genome_id
        bundle_dir.mkdir(parents=True, exist_ok=True)
        manifest = {
            "genome_id": genome_id,
            "sample_metadata": stage_file(Path(args.sample_metadata), bundle_dir, "sample_metadata.csv"),
            "reference_sheet": stage_file(Path(args.reference_sheet), bundle_dir, "reference_sheet.csv"),
        }
        for field in REQUIRED_FIELDS:
            source = genome_files[field]
            manifest[field] = stage_file(source, bundle_dir, source.name)
        write_manifest(bundle_dir / "dashboard_input_manifest.tsv", manifest)


def get_args():
    parser = argparse.ArgumentParser(description="Build per-genome LoRA dashboard input manifests")
    parser.add_argument("--run_id", required=True)
    parser.add_argument("--sample_metadata", required=True)
    parser.add_argument("--reference_sheet", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--sample_counts", nargs="*", default=[])
    parser.add_argument("--input_read_counts", nargs="*", default=[])
    parser.add_argument("--model_counts", nargs="*", default=[])
    parser.add_argument("--sqanti3_counts", nargs="*", default=[])
    parser.add_argument("--full_length_stats", nargs="*", default=[])
    parser.add_argument("--data_checks", nargs="*", default=[])
    return parser.parse_args()


def add_simple_files(by_genome, run_id, paths, field, suffix):
    prefix = f"{run_id}_"
    for raw_path in paths:
        path = checked_path(raw_path)
        name = path.name
        if not name.startswith(prefix) or not name.endswith(suffix):
            raise SystemExit(f"Cannot infer genome from {field} filename: {path}")
        genome_id = name[len(prefix):-len(suffix)]
        add_metric(by_genome, genome_id, field, path)


def add_stage_files(by_genome, run_id, paths, metric_suffix):
    prefix = f"{run_id}_"
    for raw_path in paths:
        path = checked_path(raw_path)
        name = path.name
        matched = False
        for stage in ("prefiltered", "postfiltered"):
            suffix = f"_{stage}_{metric_suffix}"
            if name.startswith(prefix) and name.endswith(suffix):
                genome_id = name[len(prefix):-len(suffix)]
                add_metric(by_genome, genome_id, f"{stage}_{metric_suffix.replace('.tsv', '')}", path)
                matched = True
                break
        if not matched:
            raise SystemExit(f"Cannot infer genome/stage from dashboard metric filename: {path}")


def add_metric(by_genome, genome_id, field, path):
    genome_files = by_genome.setdefault(genome_id, {})
    if field in genome_files:
        raise SystemExit(f"Duplicate dashboard metric for {genome_id}: {field}")
    genome_files[field] = path


def checked_path(raw_path):
    path = Path(raw_path)
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"Dashboard metric file is missing or empty: {path}")
    return path


def stage_file(source, bundle_dir, name):
    source = source.resolve()
    destination = bundle_dir / name
    if destination.exists() or destination.is_symlink():
        destination.unlink()
    try:
        os.symlink(source, destination)
    except OSError:
        shutil.copy2(source, destination)
    return destination.name


def write_manifest(path, manifest):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(manifest)


if __name__ == "__main__":
    main()
