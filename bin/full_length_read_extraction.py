#!/usr/bin/env python
import argparse
import csv
import os
import re
import subprocess
import tempfile
from collections import defaultdict


def main():
    args = get_args()
    model_multi = load_model_chains(args.models_intron_chains)
    model_single = load_model_chains(args.models_single_exons)
    read_multi = load_read_chains(args.reads_intron_chains)
    read_single = load_read_chains(args.reads_single_exons)

    full_length_reads = {}
    full_length_seen = defaultdict(set)
    full_length_multi_model_reads = defaultdict(set)
    full_length_single_model_reads = defaultdict(set)
    stats_paths = open_stat_tempfiles()

    with open(args.all_contains, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t", 2)
            if len(fields) < 3:
                continue
            chromosome, model, reads = fields
            for read in iter_comma_values(reads):
                if not read:
                    continue
                stats_paths["total_reads"].write(read + "\n")
                if read in read_multi:
                    stats_paths["total_multi"].write(read + "\n")
                    if read_multi[read] in model_multi.get(model, set()):
                        add_full_length(full_length_reads, full_length_seen, chromosome, model, read)
                        full_length_multi_model_reads[model].add(read)
                        stats_paths["fl_multi"].write(read + "\n")
                else:
                    stats_paths["total_single"].write(read + "\n")
                    chain = read_single.get(read)
                    if chain and chain in model_single.get(model, set()):
                        add_full_length(full_length_reads, full_length_seen, chromosome, model, read)
                        full_length_single_model_reads[model].add(read)
                        stats_paths["fl_single"].write(read + "\n")

    close_tempfiles(stats_paths)
    write_stats(args.output_prefix, args.all_contains, stats_paths, full_length_multi_model_reads, full_length_single_model_reads)
    write_full_length_contains(args.output_prefix, full_length_reads)
    remove_tempfiles(stats_paths)


def get_args():
    parser = argparse.ArgumentParser(description="Extract full-length read support from LoRA intron-chain tables")
    parser.add_argument("models_intron_chains")
    parser.add_argument("models_single_exons")
    parser.add_argument("reads_intron_chains")
    parser.add_argument("reads_single_exons")
    parser.add_argument("all_contains")
    parser.add_argument("output_prefix")
    return parser.parse_args()


def load_model_chains(path):
    chains = defaultdict(set)
    with open(path, encoding="utf-8") as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if len(fields) >= 3:
                chains[fields[1]].add(fields[2])
    return chains


def load_read_chains(path):
    chains = {}
    with open(path, encoding="utf-8") as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if len(fields) >= 3:
                chains[fields[1]] = fields[2]
    return chains


def iter_comma_values(text):
    start = 0
    while True:
        index = text.find(",", start)
        if index == -1:
            yield text[start:]
            return
        yield text[start:index]
        start = index + 1


def add_full_length(full_length_reads, full_length_seen, chromosome, model, read):
    if read in full_length_seen[model]:
        return
    full_length_seen[model].add(read)
    entry = full_length_reads.setdefault(model, {"chromosome": chromosome, "reads": []})
    entry["reads"].append(read)


def open_stat_tempfiles():
    paths = {}
    for name in ("total_reads", "total_multi", "total_single", "fl_multi", "fl_single"):
        handle = tempfile.NamedTemporaryFile("w", prefix=f"lora_{name}.", suffix=".txt", delete=False)
        paths[name] = handle
    return paths


def close_tempfiles(handles):
    for handle in handles.values():
        handle.close()


def remove_tempfiles(handles):
    for handle in handles.values():
        try:
            os.remove(handle.name)
        except FileNotFoundError:
            pass


def unique_line_count(path):
    sort = subprocess.Popen(["sort", "-u", path], stdout=subprocess.PIPE, text=True)
    wc = subprocess.run(["wc", "-l"], stdin=sort.stdout, stdout=subprocess.PIPE, text=True, check=True)
    if sort.stdout:
        sort.stdout.close()
    sort.wait()
    if sort.returncode != 0:
        raise RuntimeError(f"sort -u failed for {path}")
    return int(wc.stdout.strip().split()[0]) if wc.stdout.strip() else 0


def count_models(path):
    models = set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t", 2)
            if len(fields) >= 2:
                models.add(fields[1])
    return len(models)


def write_stats(output_prefix, all_contains, tempfiles, full_length_multi_model_reads, full_length_single_model_reads):
    total_reads = unique_line_count(tempfiles["total_reads"].name)
    total_multi = unique_line_count(tempfiles["total_multi"].name)
    total_single = unique_line_count(tempfiles["total_single"].name)
    fl_multi = unique_line_count(tempfiles["fl_multi"].name)
    fl_single = unique_line_count(tempfiles["fl_single"].name)
    rows = [[
        "Sample",
        "Total Reads",
        "Total Multi-Exon Reads",
        "Total Single-Exon Reads",
        "Total Full Length Reads",
        "Total Full Length Multi-Exon Reads",
        "Total Full Length Single-Exon Reads",
        "Total Transcripts",
        "Total Multi-Exon Transcripts with more that 1 FL read",
        "Total Single-Exon Transcripts with more that 1 FL read",
    ], [
        output_prefix,
        total_reads,
        total_multi,
        total_single,
        fl_multi + fl_single,
        fl_multi,
        fl_single,
        count_models(all_contains),
        sum(1 for reads in full_length_multi_model_reads.values() if len(reads) > 1),
        sum(1 for reads in full_length_single_model_reads.values() if len(reads) > 1),
    ]]
    with open(f"{output_prefix}_stats.tsv", "w", newline="", encoding="utf-8") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)


def write_full_length_contains(output_prefix, full_length_reads):
    rows = []
    for model, entry in full_length_reads.items():
        rows.append((entry["chromosome"], model, ",".join(entry["reads"])))
    rows.sort(key=lambda row: (row[0], model_number(row[1])))
    with open(f"{output_prefix}_full_length_contains.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def model_number(model):
    match = re.search(r"_(\d+)$", model)
    return int(match.group(1)) if match else 0


if __name__ == "__main__":
    main()
