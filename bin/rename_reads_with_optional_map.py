#!/usr/bin/env python3
import argparse
import csv
import gzip
import sys


COMPACT_MAP_COLUMNS = {"original_read_id", "new_read_id"}
LEGACY_MAP_COLUMNS = {
    "source_file",
    "original_read_id",
    "original_header",
    "new_read_id",
    "record_number",
}


def main():
    args = get_args()
    supplied_map = load_supplied_map(args.provided_map)
    if args.format == "fastq":
        read_count, total_bases = rename_fastq(args, supplied_map)
    else:
        read_count, total_bases = rename_fasta(args, supplied_map)
    average_read_length = validate_average_read_length(args, read_count, total_bases)
    write_summary(args, read_count, average_read_length)


def get_args():
    parser = argparse.ArgumentParser(description="Rename LoRA reads, optionally replaying a previous read ID map")
    parser.add_argument("--input", required=True)
    parser.add_argument("--format", required=True, choices=["fastq", "fasta"])
    parser.add_argument("--accession_id", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--output_map", required=True)
    parser.add_argument("--source_file", default="")
    parser.add_argument("--summary", required=True)
    parser.add_argument("--seq_source", default="")
    parser.add_argument("--genome_id", default="")
    parser.add_argument("--tissue", default="")
    parser.add_argument("--cell_line", default="")
    parser.add_argument("--provided_map", default="NA")
    parser.add_argument("--min_average_read_length", type=int, default=200)
    return parser.parse_args()


def is_provided(value):
    return str(value or "").strip().lower() not in ["", "na", "none", "null", "."]


def open_maybe_gzip(path, mode="rt"):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def load_supplied_map(path):
    if not is_provided(path):
        return None

    by_header = {}
    by_id = {}
    new_ids = set()
    with open_maybe_gzip(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail(f"Supplied read ID map {path} is empty")
        fieldnames = set(reader.fieldnames)
        if LEGACY_MAP_COLUMNS.issubset(fieldnames):
            has_header = True
            missing = set()
        else:
            has_header = False
            missing = COMPACT_MAP_COLUMNS.difference(fieldnames)
        if missing:
            fail(f"Supplied read ID map {path} is missing required columns: {', '.join(sorted(missing))}")

        for line_number, row in enumerate(reader, start=2):
            original_header = row.get("original_header", "").strip() if has_header else ""
            original_read_id = row.get("original_read_id", "").strip()
            new_read_id = row.get("new_read_id", "").strip()

            if not original_read_id or not new_read_id:
                fail(f"Supplied read ID map {path} has a blank original/new read ID on line {line_number}")
            if original_header and original_header in by_header:
                fail(f"Supplied read ID map {path} has duplicate original_header '{original_header}'")
            if original_read_id in by_id:
                fail(f"Supplied read ID map {path} has duplicate original_read_id '{original_read_id}'")
            if new_read_id in new_ids:
                fail(f"Supplied read ID map {path} has duplicate new_read_id '{new_read_id}'")

            if original_header:
                by_header[original_header] = row
            by_id[original_read_id] = row
            new_ids.add(new_read_id)

    return {"path": path, "by_header": by_header, "by_id": by_id}


def rename_fastq(args, supplied_map):
    record_number = 0
    total_bases = 0
    with open_maybe_gzip(args.input, "rt") as source, open(args.output, "w") as out, open_map(args.output_map) as map_out:
        while True:
            header_line = source.readline()
            if header_line == "":
                break
            sequence = source.readline()
            plus = source.readline()
            quality = source.readline()
            if not sequence or not plus or not quality:
                fail(f"FASTQ input {args.input} ended part-way through a record")
            if not header_line.startswith("@"):
                fail(f"FASTQ input {args.input} has a record header that does not start with @ at record {record_number + 1}")
            if not plus.startswith("+"):
                fail(f"FASTQ input {args.input} has a plus line that does not start with + at record {record_number + 1}")

            record_number += 1
            header = header_line[1:].rstrip("\n\r")
            original_id, new_id = resolve_new_id(args, supplied_map, header, record_number)
            map_out.writerow([original_id, new_id])
            out.write(f">{new_id}\n")
            out.write(sequence)
            total_bases += len(sequence.rstrip("\n\r"))
    if record_number == 0:
        fail(f"FASTQ input {args.input} did not contain any records")
    return record_number, total_bases


def rename_fasta(args, supplied_map):
    record_number = 0
    total_bases = 0
    with open_maybe_gzip(args.input, "rt") as source, open(args.output, "w") as out, open_map(args.output_map) as map_out:
        for line in source:
            if line.startswith(">"):
                record_number += 1
                header = line[1:].rstrip("\n\r")
                original_id, new_id = resolve_new_id(args, supplied_map, header, record_number)
                map_out.writerow([original_id, new_id])
                out.write(f">{new_id}\n")
            else:
                out.write(line)
                total_bases += len(line.rstrip("\n\r"))
    if record_number == 0:
        fail(f"FASTA input {args.input} did not contain any records")
    return record_number, total_bases


def resolve_new_id(args, supplied_map, header, record_number):
    original_id = header.split()[0] if header.split() else header
    if supplied_map is None:
        return original_id, f"{args.accession_id}_{record_number}"

    row = supplied_map["by_header"].get(header)
    if row is None:
        row = supplied_map["by_id"].get(original_id)
    if row is None:
        fail(f"Read '{original_id}' from {args.input} is missing from supplied read ID map {supplied_map['path']}")
    return original_id, row["new_read_id"].strip()


def source_file(args):
    return args.source_file if is_provided(args.source_file) else args.input


def metadata_value(value):
    return str(value).strip() if is_provided(value) else "not_provided"


def validate_average_read_length(args, read_count, total_bases):
    average_read_length = total_bases / read_count
    if average_read_length <= args.min_average_read_length:
        fail(
            f"{args.format.upper()} input {args.input} has average read length "
            f"{average_read_length:.0f}, which is not greater than {args.min_average_read_length}"
        )
    return round(average_read_length, 2)


def open_map(path):
    handle = open(path, "w", newline="")
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["original_read_id", "new_read_id"])

    class MapWriter:
        def __enter__(self):
            return writer

        def __exit__(self, exc_type, exc_value, traceback):
            handle.close()

    return MapWriter()


def write_summary(args, read_count, average_read_length):
    with open(args.summary, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["source_file", "format", "seq_source", "genome_id", "tissue", "cell_line", "read_count", "average_read_length"])
        writer.writerow([
            source_file(args),
            args.format,
            metadata_value(args.seq_source),
            metadata_value(args.genome_id),
            metadata_value(args.tissue),
            metadata_value(args.cell_line),
            read_count,
            average_read_length,
        ])


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
