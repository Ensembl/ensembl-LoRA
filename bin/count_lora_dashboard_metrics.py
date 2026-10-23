#!/usr/bin/env python
import argparse
import csv
import json
import os
import re
import subprocess
import tempfile
from collections import Counter, defaultdict


def main():
    args = get_args()
    os.makedirs(args.output_dir, exist_ok=True)
    if args.command == "samples":
        count_samples(args)
    elif args.command == "input_reads":
        count_input_reads(args)
    elif args.command == "models":
        count_models(args)
    elif args.command == "sqanti3":
        count_sqanti3(args)


def get_args():
    parser = argparse.ArgumentParser(description="Pre-calculate reusable LoRA dashboard metrics")
    subparsers = parser.add_subparsers(dest="command", required=True)

    samples = subparsers.add_parser("samples")
    samples.add_argument("--run_id", required=True)
    samples.add_argument("--genome_id", required=True)
    samples.add_argument("--sample_metadata", required=True)
    samples.add_argument("--output_dir", required=True)

    input_reads = subparsers.add_parser("input_reads")
    input_reads.add_argument("--run_id", required=True)
    input_reads.add_argument("--genome_id", required=True)
    input_reads.add_argument("--read_summary", required=True)
    input_reads.add_argument("--output_dir", required=True)

    models = subparsers.add_parser("models")
    models.add_argument("--run_id", required=True)
    models.add_argument("--genome_id", required=True)
    models.add_argument("--stage", required=True, choices=["prefiltered", "postfiltered"])
    models.add_argument("--gtf", required=True)
    models.add_argument("--contains", required=True)
    models.add_argument("--output_dir", required=True)

    sqanti3 = subparsers.add_parser("sqanti3")
    sqanti3.add_argument("--run_id", required=True)
    sqanti3.add_argument("--genome_id", required=True)
    sqanti3.add_argument("--stage", required=True, choices=["prefiltered", "postfiltered"])
    sqanti3.add_argument("--sqanti3_dir", required=True)
    sqanti3.add_argument("--run_sqanti3", required=True)
    sqanti3.add_argument("--output_dir", required=True)

    return parser.parse_args()


def count_samples(args):
    rows = read_dicts(args.sample_metadata, ",")
    rows = [row for row in rows if row.get("genome_id") == args.genome_id]
    samples = grouped_samples(rows)
    method_counts = Counter(value_or_unknown(first_provided(sample_rows, "seq_source")) for sample_rows in samples.values())
    tissue_counts = Counter(clean(first_provided(sample_rows, "tissue")) for sample_rows in samples.values() if provided(first_provided(sample_rows, "tissue")))
    cell_line_counts = Counter(clean(first_provided(sample_rows, "cell_line")) for sample_rows in samples.values() if provided(first_provided(sample_rows, "cell_line")))
    missing = sum(1 for sample_rows in samples.values() if not provided(first_provided(sample_rows, "tissue")) and not provided(first_provided(sample_rows, "cell_line")))
    grouped_file_count = sum(len(sample_rows) for sample_rows in samples.values() if len(sample_rows) > 1)
    output = os.path.join(args.output_dir, f"{args.run_id}_{args.genome_id}_sample_counts.tsv")
    write_rows(output, [
        ["metric_type", "group", "label", "count"],
        ["sample", "", "input_files", len(rows)],
        ["sample", "", "total_samples", len(samples)],
        ["sample", "", "unique_tissues", len(tissue_counts)],
        ["sample", "", "unique_cell_lines", len(cell_line_counts)],
        ["sample", "", "samples_missing_tissue_cell_line", missing],
        ["sample_grouping", "", "grouped_input_files", grouped_file_count],
        ["sample_grouping", "", "collapsed_sample_groups", sum(1 for sample_rows in samples.values() if len(sample_rows) > 1)],
        *[["method", "", label, count] for label, count in method_counts.most_common()],
        *[["tissue", "", label, count] for label, count in tissue_counts.most_common()],
        *[["cell_line", "", label, count] for label, count in cell_line_counts.most_common()],
    ])


def grouped_samples(rows):
    grouped = defaultdict(list)
    for row in rows:
        sample_id = row.get("sample_id")
        if not provided(sample_id):
            sample_id = row.get("reads", "")
        grouped[sample_id].append(row)
    return grouped


def first_provided(rows, column):
    for row in rows:
        value = row.get(column, "")
        if provided(value):
            return value
    return ""


def count_input_reads(args):
    total = 0
    total_bases = 0.0
    average_lengths = []
    method_counts = Counter()
    tissue_counts = Counter()
    cell_line_counts = Counter()
    tissue_cell_counts = Counter()
    missing_tissue = 0
    missing_cell_line = 0
    missing_sample_metadata = 0

    for row in iter_dicts(args.read_summary, "\t"):
        if row.get("genome_id") != args.genome_id:
            continue
        read_count = parse_int(row.get("read_count"), 0)
        average_read_length = parse_float(row.get("average_read_length"), 0.0)
        total += read_count
        if read_count > 0 and average_read_length > 0:
            total_bases += read_count * average_read_length
            average_lengths.append(average_read_length)
        source_file = row.get("source_file", "")
        method = value_or_unknown(row.get("seq_source"))
        tissue = metadata_value(row.get("tissue"))
        cell_line = metadata_value(row.get("cell_line"))

        if tissue == "not_provided":
            missing_tissue += read_count
        if cell_line == "not_provided":
            missing_cell_line += read_count
        if not provided(source_file):
            missing_sample_metadata += read_count

        method_counts[method] += read_count
        tissue_counts[tissue] += read_count
        cell_line_counts[cell_line] += read_count
        tissue_cell_counts[f"{tissue}|{cell_line}"] += read_count

    weighted_average = round(total_bases / total, 2) if total and total_bases else 0
    output = os.path.join(args.output_dir, f"{args.run_id}_{args.genome_id}_input_read_counts.tsv")
    rows = [
        ["metric_type", "group", "label", "value", "source"],
        ["total", "genome", "input_reads", total, "read_summary"],
        ["read_prep", "average_read_length", "weighted_mean", weighted_average, "read_summary"],
        ["read_prep", "average_read_length", "min", min(average_lengths) if average_lengths else 0, "read_summary"],
        ["read_prep", "average_read_length", "max", max(average_lengths) if average_lengths else 0, "read_summary"],
        ["metadata", "missing_tissue", "input_reads", missing_tissue, "read_summary"],
        ["metadata", "missing_cell_line", "input_reads", missing_cell_line, "read_summary"],
        ["metadata", "missing_sample_metadata", "input_reads", missing_sample_metadata, "read_summary"],
        *[["method", label, "input_reads", count, "read_summary"] for label, count in method_counts.most_common()],
        *[["tissue", label, "input_reads", count, "read_summary"] for label, count in tissue_counts.most_common()],
        *[["cell_line", label, "input_reads", count, "read_summary"] for label, count in cell_line_counts.most_common()],
        *[["tissue_cell", label, "input_reads", count, "read_summary"] for label, count in tissue_cell_counts.most_common()],
    ]
    write_rows(output, rows)


def count_models(args):
    models = parse_gtf(args.gtf)
    lengths = [model["mature_length"] for model in models.values()]
    single_exon = sum(1 for model in models.values() if model["exons"] == 1)
    support_counts, unique_support_reads = parse_contains_counts(args.contains, set(models), args.output_dir)
    with_support = sum(1 for count in support_counts.values() if count > 0)
    support_bins = Counter(support_bin(support_counts.get(transcript_id, 0)) for transcript_id in models)
    output = os.path.join(args.output_dir, f"{args.run_id}_{args.genome_id}_{args.stage}_model_counts.tsv")
    rows = [
        ["metric_type", "group", "label", "value"],
        ["model", args.stage, "transcript_count", len(models)],
        ["model", args.stage, "single_exon", single_exon],
        ["model", args.stage, "multi_exon", len(models) - single_exon],
        ["model", args.stage, "median_length", median(lengths)],
        ["model", args.stage, "min_length", min(lengths) if lengths else 0],
        ["model", args.stage, "max_length", max(lengths) if lengths else 0],
        ["model", args.stage, "with_support", with_support],
        ["model", args.stage, "unique_support_reads", unique_support_reads],
        *[["support_bin", args.stage, label, count] for label, count in sorted(support_bins.items())],
    ]
    write_rows(output, rows)


def count_sqanti3(args):
    prefix = f"{args.run_id}_{args.genome_id}_{args.stage}"
    classification = find_file(args.sqanti3_dir, f"{prefix}_classification.txt")
    output = os.path.join(args.output_dir, f"{args.run_id}_{args.genome_id}_{args.stage}_sqanti3_counts.tsv")
    rows = [["stage", "summary_type", "group", "label", "count"]]
    if not as_bool(args.run_sqanti3) or not classification:
        rows.append([args.stage, "available", "", "false", 0])
        write_rows(output, rows)
        return

    records, column_counts, non_canonical, subcategory_counts = count_sqanti3_classification(classification)
    rows.append([args.stage, "available", "", "true", records])
    for column in ("structural_category", "all_canonical", "coding", "CDS_type", "predicted_NMD", "RTS_stage", "within_CAGE_peak", "within_polyA_site"):
        for label, count in column_counts[column].most_common():
            rows.append([args.stage, column, "", label, count])
    rows.append([args.stage, "non_canonical", "", "non_canonical", non_canonical])
    for (structural_category, subcategory), count in sorted(subcategory_counts.items()):
        rows.append([args.stage, "subcategory", structural_category, subcategory, count])
    write_rows(output, rows)


def parse_gtf(path):
    transcripts = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            attrs = parse_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue
            start = int(fields[3])
            end = int(fields[4])
            entry = transcripts.setdefault(transcript_id, {"exons": 0, "mature_length": 0})
            entry["exons"] += 1
            entry["mature_length"] += end - start + 1
    return transcripts


def parse_attributes(text):
    return {key: value for key, value in re.findall(r'(\S+)\s+"([^"]*)"', text)}


def parse_contains_counts(path, transcript_ids, output_dir):
    support_counts = {}
    if not os.path.exists(path):
        return support_counts, 0

    unique_count = 0
    with tempfile.NamedTemporaryFile("w", dir=output_dir, prefix="support_reads.", suffix=".txt", delete=False) as read_handle:
        read_ids_path = read_handle.name
        with open(path) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t", 2)
                if len(fields) < 3:
                    continue
                transcript_id = fields[1]
                if transcript_id not in transcript_ids:
                    continue
                read_count = 0
                for read_id in iter_comma_values(fields[2]):
                    if not read_id:
                        continue
                    read_count += 1
                    read_handle.write(f"{read_id}\n")
                support_counts[transcript_id] = read_count
    try:
        unique_count = count_unique_lines(read_ids_path)
    finally:
        try:
            os.remove(read_ids_path)
        except FileNotFoundError:
            pass
    return support_counts, unique_count


def iter_comma_values(text):
    start = 0
    while True:
        index = text.find(",", start)
        if index == -1:
            yield text[start:]
            return
        yield text[start:index]
        start = index + 1


def count_unique_lines(path):
    if not os.path.exists(path):
        return 0
    sort = subprocess.Popen(["sort", "-u", path], stdout=subprocess.PIPE, text=True)
    wc = subprocess.run(["wc", "-l"], stdin=sort.stdout, stdout=subprocess.PIPE, text=True, check=True)
    if sort.stdout:
        sort.stdout.close()
    sort.wait()
    if sort.returncode != 0:
        raise RuntimeError(f"sort -u failed for temporary read ID file {path}")
    return int(wc.stdout.strip().split()[0]) if wc.stdout.strip() else 0


def count_sqanti3_classification(path):
    columns = ("structural_category", "all_canonical", "coding", "CDS_type", "predicted_NMD", "RTS_stage", "within_CAGE_peak", "within_polyA_site")
    column_counts = {column: Counter() for column in columns}
    subcategory_counts = Counter()
    non_canonical = 0
    record_count = 0
    canonical_values = {"canonical", "true", "t", "yes", "1"}
    for row in iter_dicts(path, "\t"):
        record_count += 1
        for column in columns:
            if column in row:
                column_counts[column][row.get(column, "NA") or "NA"] += 1
        if "all_canonical" in row and (row.get("all_canonical", "") or "").strip().lower() not in canonical_values:
            non_canonical += 1
        if "structural_category" in row and "subcategory" in row:
            subcategory_counts[(row.get("structural_category", "NA") or "NA", row.get("subcategory", "NA") or "NA")] += 1
    return record_count, column_counts, non_canonical, subcategory_counts


def find_file(directory, filename):
    if not os.path.isdir(directory):
        return ""
    for root, _dirs, files in os.walk(directory):
        if filename in files:
            return os.path.join(root, filename)
    return ""


def read_dicts(path, delimiter):
    if not path or not os.path.exists(path):
        return []
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def iter_dicts(path, delimiter):
    if not path or not os.path.exists(path):
        return
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            yield row


def count_column(rows, column):
    if not rows or column not in rows[0]:
        return Counter()
    return Counter(row.get(column, "NA") or "NA" for row in rows)


def count_grouped(rows, group_column, value_column):
    grouped = defaultdict(Counter)
    if not rows or group_column not in rows[0] or value_column not in rows[0]:
        return []
    for row in rows:
        grouped[row.get(group_column, "NA") or "NA"][row.get(value_column, "NA") or "NA"] += 1
    output = []
    for group, counts in sorted(grouped.items()):
        for label, count in counts.most_common():
            output.append({"structural_category": group, "subcategory": label, "count": count})
    return output


def count_non_canonical(rows):
    if not rows or "all_canonical" not in rows[0]:
        return 0
    canonical_values = {"canonical", "true", "t", "yes", "1"}
    return sum(1 for row in rows if (row.get("all_canonical", "") or "").strip().lower() not in canonical_values)


def support_bin(count):
    if count <= 0:
        return "0"
    if count == 1:
        return "1"
    if count <= 5:
        return "2-5"
    if count <= 20:
        return "6-20"
    return ">20"


def median(values):
    if not values:
        return 0
    values = sorted(values)
    mid = len(values) // 2
    if len(values) % 2:
        return values[mid]
    return round((values[mid - 1] + values[mid]) / 2, 2)


def provided(value):
    return clean(value).lower() not in {"", "na", "none", "null", "."}


def clean(value):
    return "" if value is None else str(value).strip()


def value_or_unknown(value):
    return clean(value) if provided(value) else "unknown"


def metadata_value(value):
    return clean(value) if provided(value) else "not_provided"


def parse_int(value, default=0):
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return default


def parse_float(value, default=0.0):
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return default


def unique_by_basename(rows, column):
    by_basename = defaultdict(list)
    for row in rows:
        value = row.get(column, "")
        if value:
            by_basename[os.path.basename(value)].append(row)
    return {
        basename: matches[0]
        for basename, matches in by_basename.items()
        if len(matches) == 1
    }


def unique_by_legacy_basename(rows, column):
    by_basename = defaultdict(list)
    for row in rows:
        value = row.get(column, "")
        if value:
            by_basename[legacy_lora_intermediate_basename(value)].append(row)
    return {
        basename: matches[0]
        for basename, matches in by_basename.items()
        if basename and len(matches) == 1
    }


def legacy_lora_intermediate_basename(path):
    basename = os.path.basename(path)
    basename = basename.replace(".fastq.gz", "_fastq.gz")
    basename = basename.replace(".fasta.gz", "_fasta.gz")
    return basename


def normalize_lora_intermediate_basename(path):
    basename = os.path.basename(path)
    basename = basename.replace("_fastq.gz", ".fastq.gz")
    basename = basename.replace("_fasta.gz", ".fasta.gz")
    return basename


def as_bool(value):
    return str(value).strip().lower() in {"true", "1", "yes"}


def write_rows(path, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerows(rows)


def write_dict_rows(path, rows, fieldnames):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


if __name__ == "__main__":
    main()
