#!/usr/bin/env python3
import argparse
import csv
import os
from collections import defaultdict


CHECK_HEADER = ["stage", "check_name", "status", "observed", "expected", "message", "source_file"]


def main():
    args = get_args()
    output_dir = os.path.join(args.output_dir, "data_checks")
    os.makedirs(output_dir, exist_ok=True)
    for genome_id in genome_ids(args):
        checks = []
        checks.extend(input_checks(args, genome_id))
        checks.extend(read_prep_checks(args, genome_id))
        checks.extend(model_checks(args, genome_id))
        checks.extend(full_length_checks(args, genome_id))
        checks.extend(sqanti3_checks(args, genome_id))
        checks.extend(bam_checks(args, genome_id))
        checks.extend(dashboard_input_checks(args))
        output = os.path.join(output_dir, f"{args.run_id}_{genome_id}_data_checks.tsv")
        write_rows(output, [CHECK_HEADER] + checks)


def get_args():
    parser = argparse.ArgumentParser(description="Build reusable LoRA data-check metrics")
    parser.add_argument("--run_id", required=True)
    parser.add_argument("--sample_metadata", required=True)
    parser.add_argument("--sample_counts", nargs="*", default=[])
    parser.add_argument("--input_read_counts", nargs="*", default=[])
    parser.add_argument("--model_counts", nargs="*", default=[])
    parser.add_argument("--sqanti3_counts", nargs="*", default=[])
    parser.add_argument("--full_length_stats", nargs="*", default=[])
    parser.add_argument("--bam_stats", nargs="*", default=[])
    parser.add_argument("--output_dir", required=True)
    return parser.parse_args()


def genome_ids(args):
    ids = {row.get("genome_id") for row in read_dicts(args.sample_metadata, ",") if provided(row.get("genome_id"))}
    for paths in (args.sample_counts, args.input_read_counts, args.model_counts, args.sqanti3_counts, args.full_length_stats):
        for path in paths:
            inferred = infer_genome_from_metric(args.run_id, path)
            if inferred:
                ids.add(inferred)
    return sorted(ids)


def input_checks(args, genome_id):
    rows = [row for row in read_dicts(args.sample_metadata, ",") if row.get("genome_id") == genome_id]
    paths = [row.get("reads", "") for row in rows]
    samples = grouped_samples(rows)
    duplicates = len(paths) - len(set(paths))
    missing_tissue = sum(1 for row in rows if not provided(row.get("tissue")))
    missing_cell_line = sum(1 for row in rows if not provided(row.get("cell_line")))
    missing_sample_id = sum(1 for row in rows if not provided(row.get("sample_id")))
    grouped_files = sum(len(sample_rows) for sample_rows in samples.values() if len(sample_rows) > 1)
    collapsed_groups = sum(1 for sample_rows in samples.values() if len(sample_rows) > 1)
    conflicts = sample_metadata_conflicts(samples)
    checks = [
        check("Inputs", "sample_rows", "passed" if rows else "failed", len(rows), ">0", "Sample metadata rows for this genome", args.sample_metadata),
        check("Inputs", "input_files", "passed" if rows else "failed", len(rows), ">0", "Physical input files for this genome", args.sample_metadata),
        check("Inputs", "sample_groups", "passed" if samples else "failed", len(samples), ">0", "Dashboard sample groups for this genome", args.sample_metadata),
        check("Inputs", "grouped_input_files", "passed", grouped_files, ">=0", "Input files represented inside multi-file sample groups", args.sample_metadata),
        check("Inputs", "collapsed_sample_groups", "passed", collapsed_groups, ">=0", "Sample groups built from more than one input file", args.sample_metadata),
        check("Inputs", "duplicate_read_paths", "warning" if duplicates else "passed", duplicates, "0", "Duplicate read paths in sample metadata", args.sample_metadata),
        check("Inputs", "missing_tissue_values", "warning" if missing_tissue else "passed", missing_tissue, "0", "Rows without tissue metadata", args.sample_metadata),
        check("Inputs", "missing_cell_line_values", "warning" if missing_cell_line else "passed", missing_cell_line, "0", "Rows without cell-line metadata", args.sample_metadata),
        check("Inputs", "missing_sample_id_values", "warning" if missing_sample_id else "passed", missing_sample_id, "0", "Rows without an explicit or inferred sample_id", args.sample_metadata),
        check("Inputs", "sample_group_metadata_conflicts", "warning" if conflicts else "passed", len(conflicts), "0", "Sample groups with conflicting method/tissue/cell-line/species/map_key/genome metadata", args.sample_metadata),
    ]
    for sample_id, columns in conflicts[:20]:
        checks.append(check("Inputs", f"sample_group_conflict:{sample_id}", "warning", ",".join(columns), "consistent metadata", "Conflicting metadata within grouped input files", args.sample_metadata))
    return checks


def grouped_samples(rows):
    grouped = defaultdict(list)
    for row in rows:
        sample_id = row.get("sample_id") if provided(row.get("sample_id")) else row.get("reads", "")
        grouped[sample_id].append(row)
    return grouped


def sample_metadata_conflicts(samples):
    checked_columns = ["seq_source", "tissue", "cell_line", "species", "map_key", "genome_id"]
    conflicts = []
    for sample_id, rows in sorted(samples.items()):
        conflicting_columns = []
        for column in checked_columns:
            values = {row.get(column, "") for row in rows if provided(row.get(column, ""))}
            if len(values) > 1:
                conflicting_columns.append(column)
        if conflicting_columns:
            conflicts.append((sample_id, conflicting_columns))
    return conflicts


def read_prep_checks(args, genome_id):
    input_counts = genome_file(args.input_read_counts, args.run_id, genome_id, "_input_read_counts.tsv")
    rows = read_dicts(input_counts, "\t")
    total = metric_value(rows, "total", "genome", "input_reads", "value")
    average_length = metric_value(rows, "read_prep", "average_read_length", "weighted_mean", "value")
    min_length = metric_value(rows, "read_prep", "average_read_length", "min", "value")
    max_length = metric_value(rows, "read_prep", "average_read_length", "max", "value")
    checks = [
        check("Read Prep", "input_read_count", "passed" if int_or_zero(total) > 0 else "failed", total, ">0", "Reads counted from read summaries", input_counts),
        check("Read Prep", "average_read_length", "passed" if float_or_zero(average_length) > 200 else "warning", average_length, ">200", "Weighted mean read length from read summaries", input_counts),
        check("Read Prep", "minimum_file_average_read_length", "passed" if float_or_zero(min_length) > 200 else "warning", min_length, ">200", "Lowest per-file average read length", input_counts),
        check("Read Prep", "maximum_file_average_read_length", "passed" if float_or_zero(max_length) > 0 else "warning", max_length, ">0", "Highest per-file average read length", input_counts),
    ]
    for row in rows:
        if row.get("metric_type") == "metadata" and int_or_zero(row.get("value")) > 0:
            checks.append(check("Read Prep", row.get("group", "metadata"), "warning", row.get("value", "0"), "0", "Read metadata missing for some reads", input_counts))
    return checks


def model_checks(args, genome_id):
    checks = []
    for stage in ("prefiltered", "postfiltered"):
        path = stage_file(args.model_counts, args.run_id, genome_id, stage, "model_counts.tsv")
        rows = read_dicts(path, "\t")
        transcript_count = label_value(rows, "model", stage, "transcript_count")
        single_exon = label_value(rows, "model", stage, "single_exon")
        multi_exon = label_value(rows, "model", stage, "multi_exon")
        status = "passed" if int_or_zero(transcript_count) > 0 else "failed"
        checks.append(check("Transcript Models", f"{stage}_transcript_count", status, transcript_count, ">0", f"{stage} transcript model count", path))
        checks.append(check("Transcript Models", f"{stage}_single_exon", "passed", single_exon, ">=0", f"{stage} single-exon transcript count", path))
        checks.append(check("Transcript Models", f"{stage}_multi_exon", "passed", multi_exon, ">=0", f"{stage} multi-exon transcript count", path))
    return checks


def full_length_checks(args, genome_id):
    path = genome_file(args.full_length_stats, args.run_id, genome_id, "_full_length_stats.tsv")
    rows = read_dicts(path, "\t")
    if not rows:
        return [check("Full-Length Filtering", "full_length_stats", "failed", "missing", "present", "Full-length stats file missing or empty", path)]
    row = rows[0]
    retained = row.get("Total Full Length Reads", "not_available")
    return [check("Full-Length Filtering", "full_length_reads", "passed" if int_or_zero(retained) > 0 else "warning", retained, ">0", "Full-length reads detected", path)]


def sqanti3_checks(args, genome_id):
    checks = []
    for stage in ("prefiltered", "postfiltered"):
        path = stage_file(args.sqanti3_counts, args.run_id, genome_id, stage, "sqanti3_counts.tsv")
        rows = read_dicts(path, "\t")
        available = any(row.get("summary_type") == "available" and row.get("label") == "true" for row in rows)
        status = "passed" if available else "skipped"
        count = next((row.get("count", "0") for row in rows if row.get("summary_type") == "available"), "0")
        checks.append(check("SQANTI3", f"{stage}_classification", status, count, "present when SQANTI3 enabled", f"{stage} SQANTI3 classification rows", path))
    return checks


def bam_checks(args, genome_id):
    checks = []
    by_stage = defaultdict(list)
    for path in args.bam_stats:
        for row in read_dicts(path, "\t"):
            if row.get("genome_id") == genome_id:
                by_stage[row.get("stage", "unknown")].append((path, row))
    for stage, rows in sorted(by_stage.items()):
        total = sum(int_or_zero(row.get("total")) for _path, row in rows)
        mapped = sum(int_or_zero(row.get("mapped")) for _path, row in rows)
        missing_index = sum(1 for path, row in rows if int_or_zero(row.get("index_size_bytes")) == 0)
        empty = sum(1 for _path, row in rows if int_or_zero(row.get("total")) == 0)
        source = rows[0][0] if rows else ""
        checks.append(check("BAM Processing", f"{stage}_bam_total_records", "passed" if total > 0 else "warning", total, ">0", f"{stage} BAM records across checked files", source))
        checks.append(check("Mapping", f"{stage}_mapped_records", "passed" if mapped > 0 else "warning", mapped, ">0", f"{stage} mapped BAM records", source))
        checks.append(check("BAM Processing", f"{stage}_missing_indexes", "failed" if missing_index else "passed", missing_index, "0", f"{stage} BAM files missing non-empty indexes", source))
        checks.append(check("BAM Processing", f"{stage}_empty_bams", "warning" if empty else "passed", empty, "0", f"{stage} empty BAM files", source))
    if not by_stage:
        checks.append(check("BAM Processing", "bam_stats", "skipped", "not_available", "available after BAM stages", "No BAM stats were provided", ""))
    return checks


def dashboard_input_checks(args):
    required = {
        "sample_counts": args.sample_counts,
        "input_read_counts": args.input_read_counts,
        "model_counts": args.model_counts,
        "sqanti3_counts": args.sqanti3_counts,
        "full_length_stats": args.full_length_stats,
    }
    return [
        check("Dashboard Inputs", name, "passed" if paths else "failed", len(paths), ">0", f"{name} files supplied to dashboard checks", ",".join(paths))
        for name, paths in required.items()
    ]


def infer_genome_from_metric(run_id, path):
    name = os.path.basename(path)
    prefix = f"{run_id}_"
    if not name.startswith(prefix):
        return ""
    suffixes = [
        "_sample_counts.tsv",
        "_input_read_counts.tsv",
        "_prefiltered_model_counts.tsv",
        "_postfiltered_model_counts.tsv",
        "_prefiltered_sqanti3_counts.tsv",
        "_postfiltered_sqanti3_counts.tsv",
        "_full_length_stats.tsv",
    ]
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[len(prefix):-len(suffix)]
    return ""


def check(stage, name, status, observed, expected, message, source):
    return [stage, name, status, observed, expected, message, source]


def read_dicts(path, delimiter):
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_rows(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)


def first_existing(paths, suffix):
    for path in paths:
        if path.endswith(suffix) and os.path.exists(path):
            return path
    return paths[0] if paths else ""


def genome_file(paths, run_id, genome_id, suffix):
    needle = f"{run_id}_{genome_id}{suffix}"
    for path in paths:
        if os.path.basename(path) == needle:
            return path
    return ""


def stage_file(paths, run_id, genome_id, stage, suffix):
    needle = f"{run_id}_{genome_id}_{stage}_{suffix}"
    for path in paths:
        if os.path.basename(path) == needle:
            return path
    return ""


def metric_value(rows, metric_type, group, label, field):
    for row in rows:
        if row.get("metric_type") == metric_type and row.get("group") == group and row.get("label") == label:
            return row.get(field, "0")
    return "0"


def label_value(rows, metric_type, group, label):
    for row in rows:
        if row.get("metric_type") == metric_type and row.get("group") == group and row.get("label") == label:
            return row.get("value", "0")
    return "0"


def provided(value):
    return str(value or "").strip().lower() not in {"", "na", "none", "null", "."}


def int_or_zero(value):
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return 0


def float_or_zero(value):
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return 0.0


if __name__ == "__main__":
    main()
