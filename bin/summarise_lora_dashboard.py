#!/usr/bin/env python
import argparse
import csv
import json
import os
from datetime import datetime, timezone


def main():
    args = get_args()
    apply_dashboard_input_manifest(args)
    os.makedirs(args.output_dir, exist_ok=True)

    sample_counts = read_table(args.sample_counts)
    input_read_counts = read_table(args.input_read_counts)
    data_checks = read_table(args.data_checks)
    pre_models = model_summary(args.prefiltered_model_counts)
    post_models = model_summary(args.postfiltered_model_counts)
    pre_sqanti = sqanti3_summary(args.prefiltered_sqanti3_counts)
    post_sqanti = sqanti3_summary(args.postfiltered_sqanti3_counts)

    run_summary = {
        "genome": args.genome_id,
        "run_id": args.run_id,
        "full_length_threshold": args.min_fl,
        "genomeSize": args.genome_size,
        "input_reads": grouped_metric_value(input_read_counts, "total", "genome", "input_reads", "value", "not_available"),
        "input_reads_source": grouped_metric_value(input_read_counts, "total", "genome", "input_reads", "source", "not_available"),
        "canonical_reads": pre_models.get("unique_support_reads", "not_available"),
        "input_file_count": sample_count(sample_counts, "sample", "input_files"),
        "input_sample_count": sample_count(sample_counts, "sample", "total_samples"),
        "input_tissue_count": sample_count(sample_counts, "sample", "unique_tissues"),
        "input_cell_line_count": sample_count(sample_counts, "sample", "unique_cell_lines"),
        "input_samples_missing_tissue_cell_line": sample_count(sample_counts, "sample", "samples_missing_tissue_cell_line"),
        "grouped_input_files": sample_count(sample_counts, "sample_grouping", "grouped_input_files"),
        "collapsed_sample_groups": sample_count(sample_counts, "sample_grouping", "collapsed_sample_groups"),
        "samples_per_method": count_dict(sample_counts, "method"),
        "tissue_counts": count_dict(sample_counts, "tissue"),
        "cell_line_counts": count_dict(sample_counts, "cell_line"),
        "read_counts_by_method": grouped_count_dict(input_read_counts, "method"),
        "read_counts_by_tissue": grouped_count_dict(input_read_counts, "tissue"),
        "read_counts_by_cell_line": grouped_count_dict(input_read_counts, "cell_line"),
        "read_counts_by_tissue_cell_line": grouped_count_dict(input_read_counts, "tissue_cell"),
        "read_counts_missing_metadata": grouped_count_dict(input_read_counts, "metadata"),
    }

    data = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "run_summary": run_summary,
        "models": {
            "prefiltered": pre_models,
            "postfiltered": post_models,
        },
        "sqanti3": {
            "prefiltered": pre_sqanti,
            "postfiltered": post_sqanti,
            "delta_structural_category": delta_rows(pre_sqanti.get("structural_category", []), post_sqanti.get("structural_category", [])),
            "delta_coding": delta_rows(pre_sqanti.get("coding", []), post_sqanti.get("coding", [])),
        },
        "files": {
            "sample_metadata": args.sample_metadata,
            "reference_sheet": args.reference_sheet,
            "full_length_stats": args.full_length_stats,
            "data_checks": args.data_checks,
        },
        "data_checks": data_check_summary(data_checks),
        "summary_files": {
            "sample_counts": args.sample_counts,
            "input_read_counts": args.input_read_counts,
            "prefiltered_model_counts": args.prefiltered_model_counts,
            "postfiltered_model_counts": args.postfiltered_model_counts,
            "prefiltered_sqanti3_counts": args.prefiltered_sqanti3_counts,
            "postfiltered_sqanti3_counts": args.postfiltered_sqanti3_counts,
            "data_checks": args.data_checks,
        },
    }

    run_summary_path = os.path.join(args.output_dir, f"{args.run_id}_{args.genome_id}_run_summary.json")
    dashboard_path = os.path.join(args.output_dir, f"{args.run_id}_{args.genome_id}_dashboard_data.json")
    with open(run_summary_path, "w") as handle:
        json.dump(run_summary, handle, indent=2)
    with open(dashboard_path, "w") as handle:
        json.dump(data, handle, indent=2)


def get_args():
    parser = argparse.ArgumentParser(description="Assemble precomputed LoRA dashboard metrics")
    parser.add_argument("--run_id", required=True)
    parser.add_argument("--genome_id", required=True)
    parser.add_argument("--dashboard_input_dir")
    parser.add_argument("--reference_sheet")
    parser.add_argument("--sample_metadata")
    parser.add_argument("--sample_counts")
    parser.add_argument("--input_read_counts")
    parser.add_argument("--prefiltered_model_counts")
    parser.add_argument("--postfiltered_model_counts")
    parser.add_argument("--prefiltered_sqanti3_counts")
    parser.add_argument("--postfiltered_sqanti3_counts")
    parser.add_argument("--full_length_stats")
    parser.add_argument("--data_checks")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--min_fl", required=True)
    parser.add_argument("--genome_size", required=True)
    return parser.parse_args()


def apply_dashboard_input_manifest(args):
    if not args.dashboard_input_dir:
        required = [
            "reference_sheet",
            "sample_metadata",
            "sample_counts",
            "input_read_counts",
            "prefiltered_model_counts",
            "postfiltered_model_counts",
            "prefiltered_sqanti3_counts",
            "postfiltered_sqanti3_counts",
            "full_length_stats",
            "data_checks",
        ]
        missing = [name for name in required if not getattr(args, name)]
        if missing:
            raise SystemExit(f"Missing required dashboard input arguments: {', '.join(missing)}")
        return

    input_dir = os.path.abspath(args.dashboard_input_dir)
    manifest_path = os.path.join(input_dir, "dashboard_input_manifest.tsv")
    rows = read_table(manifest_path)
    if len(rows) != 1:
        raise SystemExit(f"Expected exactly one row in dashboard input manifest: {manifest_path}")

    row = rows[0]
    for field in [
        "reference_sheet",
        "sample_metadata",
        "sample_counts",
        "input_read_counts",
        "prefiltered_model_counts",
        "postfiltered_model_counts",
        "prefiltered_sqanti3_counts",
        "postfiltered_sqanti3_counts",
        "full_length_stats",
        "data_checks",
    ]:
        value = row.get(field, "")
        if not value:
            raise SystemExit(f"Dashboard input manifest is missing field '{field}': {manifest_path}")
        path = value if os.path.isabs(value) else os.path.join(input_dir, value)
        setattr(args, field, path)


def read_table(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def model_summary(path):
    rows = read_table(path)
    summary = {
        "transcript_count": 0,
        "single_exon": 0,
        "multi_exon": 0,
        "median_length": 0,
        "min_length": 0,
        "max_length": 0,
        "with_support": 0,
        "unique_support_reads": "not_available",
        "support_bins": {},
    }
    for row in rows:
        if row.get("metric_type") == "model":
            summary[row["label"]] = number(row.get("value", 0))
        elif row.get("metric_type") == "support_bin":
            summary["support_bins"][row["label"]] = number(row.get("value", 0))
    return summary


def sqanti3_summary(path):
    rows = read_table(path)
    summary = {
        "available": False,
        "classification_rows": 0,
        "structural_category": [],
        "subcategory_by_structural_category": [],
        "all_canonical": [],
        "non_canonical": 0,
        "coding": [],
        "CDS_type": [],
        "predicted_NMD": [],
        "RTS_stage": [],
        "within_CAGE_peak": [],
        "within_polyA_site": [],
    }
    for row in rows:
        summary_type = row.get("summary_type")
        label = row.get("label", "")
        count = number(row.get("count", 0))
        if summary_type == "available":
            summary["available"] = label == "true"
            summary["classification_rows"] = count if label == "true" else 0
        elif summary_type == "subcategory":
            summary["subcategory_by_structural_category"].append({
                "structural_category": row.get("group", ""),
                "subcategory": label,
                "count": count,
            })
        elif summary_type == "non_canonical":
            summary["non_canonical"] = count
        elif summary_type in summary:
            summary[summary_type].append({"label": label, "count": count})
    return summary


def sample_count(rows, metric_type, label):
    for row in rows:
        if row.get("metric_type") == metric_type and row.get("label") == label:
            return number(row.get("count", 0))
    return 0


def count_dict(rows, metric_type):
    return {row["label"]: number(row.get("count", 0)) for row in rows if row.get("metric_type") == metric_type}


def metric_value(rows, key_field, key_value, value_field, default):
    for row in rows:
        if row.get(key_field) == key_value:
            return number(row.get(value_field, default))
    return default


def grouped_metric_value(rows, metric_type, group, label, value_field, default):
    for row in rows:
        if row.get("metric_type") == metric_type and row.get("group") == group and row.get("label") == label:
            return number(row.get(value_field, default))
    return default


def grouped_count_dict(rows, metric_type):
    return {
        row["group"]: number(row.get("value", 0))
        for row in rows
        if row.get("metric_type") == metric_type and row.get("label") == "input_reads"
    }


def delta_rows(before, after):
    before_dict = {item["label"]: item["count"] for item in before}
    after_dict = {item["label"]: item["count"] for item in after}
    return [
        {
            "label": label,
            "prefiltered": before_dict.get(label, 0),
            "postfiltered": after_dict.get(label, 0),
        }
        for label in sorted(set(before_dict) | set(after_dict))
    ]


def data_check_summary(rows):
    status_counts = {}
    stages = {}
    for row in rows:
        status = row.get("status", "unknown") or "unknown"
        stage = row.get("stage", "unknown") or "unknown"
        status_counts[status] = status_counts.get(status, 0) + 1
        stages.setdefault(stage, []).append(row)
    return {
        "status_counts": status_counts,
        "stages": stages,
        "rows": rows,
    }


def number(value):
    try:
        return int(str(value))
    except ValueError:
        try:
            return float(str(value))
        except ValueError:
            return value


if __name__ == "__main__":
    main()
