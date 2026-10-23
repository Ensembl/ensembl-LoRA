#!/usr/bin/env python3
import argparse
import csv
import os
from pathlib import Path


STAGES = ["read_prep", "mapped_bam", "bam_to_gtf", "merged_models", "full_length", "sqanti3", "dashboard"]
RESUMABLE_STAGES = ["read_prep", "mapped_bam", "merged_models", "full_length", "sqanti3", "dashboard"]
NON_GENOME_DIRS = {"dashboard", "pipeline_info", "summary", "work", "tmp", "checkpoints"}
EMPTY_VALUES = {"", "na", "none", "null", "."}


def main():
    args = get_args()
    outdir = resolve_outdir(args.outdir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sample_metadata = read_dicts(args.sample_metadata)
    references = read_reference_sheet(args.reference_sheet)
    genomes = discover_genomes(outdir)
    if not genomes:
        raise SystemExit(no_genomes_message(outdir, args.run_id))

    plan = resolve_plan(outdir, args.run_id, genomes, args.resume_from, args.resume_scope)
    write_manifests(output_dir, outdir, args.run_id, plan, sample_metadata, references)


def get_args():
    parser = argparse.ArgumentParser(description="Build LoRA published-output checkpoint loader manifests")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--run_id", required=True)
    parser.add_argument("--resume_from", required=True)
    parser.add_argument("--resume_scope", required=True, choices=["common", "per_genome"])
    parser.add_argument("--sample_metadata", required=True)
    parser.add_argument("--reference_sheet", required=True)
    parser.add_argument("--output_dir", required=True)
    return parser.parse_args()


def resolve_outdir(outdir):
    path = Path(outdir).expanduser()
    if not path.is_absolute():
        path = Path.cwd() / path
    return path.resolve()


def no_genomes_message(outdir, run_id):
    lines = [
        f"No genome output directories found under {outdir}",
        f"Scanned checkpoint outdir: {outdir}",
        f"Expected genome directories like: {outdir}/<genome_id>/",
    ]
    if not outdir.exists():
        lines.append("The checkpoint outdir does not exist.")
    elif not outdir.is_dir():
        lines.append("The checkpoint outdir exists but is not a directory.")
    else:
        entries = sorted(path.name for path in outdir.iterdir())[:30]
        lines.append("Found top-level entries: " + (", ".join(entries) if entries else "(none)"))
    outdir_text = str(outdir)
    if outdir_text.startswith("/nfs/"):
        lines.append("This looks like an EBI /nfs path. If the run used /gpfs, try the matching /gpfs path.")
    elif outdir_text.startswith("/gpfs/"):
        lines.append("This looks like an EBI /gpfs path. If the run used /nfs, check whether the same outputs exist through that mount.")
    lines.append(f"Run: lora_checkpoint_status.py --outdir {outdir} --run_id {run_id}")
    return "\n".join(lines)


def discover_genomes(outdir):
    if not outdir.exists():
        return []
    return sorted(
        path.name for path in outdir.iterdir()
        if path.is_dir() and path.name not in NON_GENOME_DIRS and not path.name.startswith(".")
    )


def resolve_plan(outdir, run_id, genomes, requested, scope):
    if requested == "bam_to_gtf":
        fail_bam_to_gtf()
    if requested == "none":
        return {}
    statuses = {genome: {stage: stage_status(outdir, run_id, genome, stage) for stage in STAGES} for genome in genomes}
    if scope == "common":
        stage = recommended_common(statuses) if requested == "auto" else requested
        validate_stage(stage, statuses, genomes)
        return {genome: stage for genome in genomes}
    plan = {}
    for genome in genomes:
        stage = recommended_for_genome(statuses[genome]) if requested == "auto" else requested
        validate_stage(stage, statuses, [genome])
        plan[genome] = stage
    return plan


def validate_stage(stage, statuses, genomes):
    if stage == "bam_to_gtf":
        fail_bam_to_gtf()
    if stage not in RESUMABLE_STAGES:
        raise SystemExit(f"Checkpoint '{stage}' is not resumable from published outputs. Run lora_checkpoint_status.py to inspect available checkpoints.")
    for genome in genomes:
        status = statuses[genome].get(stage, "missing")
        if status != "complete":
            raise SystemExit(f"Checkpoint '{stage}' is {status} for genome {genome}; cannot resume from it.")


def fail_bam_to_gtf():
    raise SystemExit("bam_to_gtf is audit-only in this LoRA version. Resume from mapped_bam to rerun BAM-to-GTF, or from merged_models if model merging completed.")


def write_manifests(output_dir, outdir, run_id, plan, sample_metadata, references):
    writers = {
        "checkpoint_plan.csv": [],
        "read_prep_fasta.csv": [],
        "mapped_bam.csv": [],
        "canonical_bam.csv": [],
        "chromosome_models.csv": [],
        "prefiltered_models.csv": [],
        "prefiltered_contains.csv": [],
        "filtered_models.csv": [],
        "filtered_chromosome_models.csv": [],
        "full_length_stats.csv": [],
        "read_summaries.csv": [],
        "sqanti3_dashboard.csv": [],
        "sample_counts.csv": [],
        "input_read_counts.csv": [],
        "prefiltered_model_counts.csv": [],
        "postfiltered_model_counts.csv": [],
        "prefiltered_sqanti3_counts.csv": [],
        "postfiltered_sqanti3_counts.csv": [],
        "data_checks.csv": [],
        "dashboard_bam_stats_inputs.csv": [],
    }
    metadata_by_genome = group_metadata_by_genome(sample_metadata)
    reference_by_genome = reference_lookup(references)
    for genome, stage in sorted(plan.items()):
        writers["checkpoint_plan.csv"].append([genome, stage])
        reference = reference_by_genome.get(genome, "")
        species = first_value(metadata_by_genome.get(genome, []), "species", "unknown")
        raw_file_count = str(len(metadata_by_genome.get(genome, [])) or 1)
        read_summary = exact_file(outdir / genome / "reads" / "read_summaries", f"{run_id}_{genome}_read_summary.tsv")
        if read_summary:
            writers["read_summaries.csv"].append([genome, read_summary])

        if stage == "read_prep":
            seq_source_by_name = seq_source_lookup(metadata_by_genome.get(genome, []))
            prepared_fastas = prepared_fasta_files(outdir / genome)
            if not prepared_fastas:
                raise SystemExit(f"Checkpoint 'read_prep' is complete in the audit only if prepared FASTA files are published; none were found for genome {genome}.")
            default_seq_source = first_value(metadata_by_genome.get(genome, []), "seq_source", "nanopore")
            for fasta in prepared_fastas:
                seq_source = seq_source_for_fasta(fasta, seq_source_by_name, default_seq_source)
                writers["read_prep_fasta.csv"].append([fasta, seq_source, reference, genome, raw_file_count, species])

        if stage == "mapped_bam":
            bam = exact_file(outdir / genome / "reads" / "all_reads" / "merged", f"{run_id}_{genome}.bam")
            bai = exact_file(outdir / genome / "reads" / "all_reads" / "merged", f"{run_id}_{genome}.bam.csi")
            writers["mapped_bam.csv"].append([bam, bai, reference, genome])

        if stage in {"merged_models", "full_length", "sqanti3", "dashboard"}:
            pre_gtf = exact_file(outdir / genome / "all_models" / "merged", f"{run_id}_{genome}.gtf")
            pre_contains = exact_file(outdir / genome / "all_models" / "merged", f"{run_id}_{genome}_contains.tsv")
            writers["prefiltered_models.csv"].append([genome, pre_gtf])
            writers["prefiltered_contains.csv"].append([genome, pre_contains])

        if stage in {"merged_models", "full_length", "sqanti3", "dashboard"}:
            for gtf, contains, chromosome_id, chromosome_count in chromosome_models(outdir / genome, run_id):
                writers["chromosome_models.csv"].append([genome, gtf, contains, chromosome_count, chromosome_id])

        if stage == "merged_models":
            for bam, bai, chromosome_id, chromosome_count in canonical_bams(outdir / genome, run_id):
                writers["canonical_bam.csv"].append([bam, bai, genome, chromosome_id, chromosome_count])

        if stage in {"full_length", "sqanti3", "dashboard"}:
            post_gtf = exact_file(outdir / genome / "full_length" / "models" / "merged", f"{run_id}_{genome}_filtered_transcripts.gtf")
            post_contains = exact_file(outdir / genome / "full_length" / "models" / "merged", f"{run_id}_{genome}_filtered_contains.tsv")
            stats = exact_file(outdir / genome / "full_length" / "stats", f"{run_id}_{genome}_full_length_stats.tsv")
            writers["filtered_models.csv"].append([genome, post_gtf, post_contains])
            writers["full_length_stats.csv"].append([genome, stats])
            for gtf, contains, chromosome_id, chromosome_count in full_length_chromosome_models(outdir / genome, run_id):
                writers["filtered_chromosome_models.csv"].append([genome, gtf, contains, chromosome_count, chromosome_id])

        if stage in {"merged_models", "full_length", "sqanti3"}:
            for row in dashboard_bam_stats_inputs(outdir / genome, run_id, genome):
                writers["dashboard_bam_stats_inputs.csv"].append(row)

        if stage in {"sqanti3", "dashboard"}:
            sqanti_status = stage_status(outdir, run_id, genome, "sqanti3")
            sqanti_dir = str(outdir / genome / "sqanti3") if sqanti_status == "complete" else str(outdir)
            writers["sqanti3_dashboard.csv"].append([genome, "true" if sqanti_status == "complete" else "false", sqanti_dir])

        if stage == "dashboard":
            metrics_dir = outdir / genome / "summary" / "metrics"
            writers["sample_counts.csv"].append([
                genome,
                exact_file(metrics_dir, f"{run_id}_{genome}_sample_counts.tsv"),
            ])
            writers["input_read_counts.csv"].append([
                genome,
                exact_file(metrics_dir, f"{run_id}_{genome}_input_read_counts.tsv"),
            ])
            writers["prefiltered_model_counts.csv"].append([
                genome,
                exact_file(metrics_dir, f"{run_id}_{genome}_prefiltered_model_counts.tsv"),
            ])
            writers["postfiltered_model_counts.csv"].append([
                genome,
                exact_file(metrics_dir, f"{run_id}_{genome}_postfiltered_model_counts.tsv"),
            ])
            writers["prefiltered_sqanti3_counts.csv"].append([
                genome,
                exact_file(metrics_dir, f"{run_id}_{genome}_prefiltered_sqanti3_counts.tsv"),
            ])
            writers["postfiltered_sqanti3_counts.csv"].append([
                genome,
                exact_file(metrics_dir, f"{run_id}_{genome}_postfiltered_sqanti3_counts.tsv"),
            ])
            checks_dir = outdir / genome / "summary" / "checks"
            writers["data_checks.csv"].append([
                genome,
                exact_file(checks_dir, f"{run_id}_{genome}_data_checks.tsv"),
            ])

    for filename, rows in writers.items():
        write_rows(output_dir / filename, rows)


def stage_status(outdir, run_id, genome, stage):
    genome_dir = outdir / genome
    if stage == "read_prep":
        read_summary = status_for_patterns(genome_dir, ["reads/read_summaries/*_read_summary.tsv"])
        prepared = "complete" if prepared_fasta_files(genome_dir) else "missing"
        if read_summary == "complete" and prepared == "complete":
            return "complete"
        if read_summary == "complete":
            return "not_resumable"
        return read_summary
    if stage == "mapped_bam":
        return status_for_patterns(genome_dir, [
            exact_or_glob(run_id, genome, "reads/all_reads/merged/{run}_{genome}.bam", "reads/all_reads/merged/*.bam"),
            exact_or_glob(run_id, genome, "reads/all_reads/merged/{run}_{genome}.bam.csi", "reads/all_reads/merged/*.bam.csi"),
        ])
    if stage == "bam_to_gtf":
        progress = status_for_patterns(genome_dir, [
            "reads/canonical/chromosomes/*.canonical.bam",
            "reads/canonical/chromosomes/*.canonical.bam.csi",
            "all_models/chromosomes/*.gtf",
        ])
        return "not_resumable" if progress in {"complete", "partial"} else progress
    if stage == "merged_models":
        return status_for_patterns(genome_dir, [
            exact_or_glob(run_id, genome, "all_models/merged/{run}_{genome}.gtf", "all_models/merged/*.gtf"),
            exact_or_glob(run_id, genome, "all_models/merged/{run}_{genome}_contains.tsv", "all_models/merged/*_contains.tsv"),
            "all_models/chromosomes/*.combined.gtf",
            "all_models/chromosomes/*.combined.contains.tsv",
            "reads/canonical/chromosomes/*.canonical.bam",
            "reads/canonical/chromosomes/*.canonical.bam.csi",
        ])
    if stage == "full_length":
        return status_for_patterns(genome_dir, [
            exact_or_glob(run_id, genome, "all_models/merged/{run}_{genome}.gtf", "all_models/merged/*.gtf"),
            exact_or_glob(run_id, genome, "all_models/merged/{run}_{genome}_contains.tsv", "all_models/merged/*_contains.tsv"),
            exact_or_glob(run_id, genome, "full_length/models/merged/{run}_{genome}_filtered_transcripts.gtf", "full_length/models/merged/*_filtered_transcripts.gtf"),
            exact_or_glob(run_id, genome, "full_length/models/merged/{run}_{genome}_filtered_contains.tsv", "full_length/models/merged/*_filtered_contains.tsv"),
            exact_or_glob(run_id, genome, "full_length/stats/{run}_{genome}_full_length_stats.tsv", "full_length/stats/*_full_length_stats.tsv"),
        ])
    if stage == "sqanti3":
        sqanti_dir = genome_dir / "sqanti3"
        if not sqanti_dir.exists():
            return "skipped"
        return status_for_patterns(genome_dir, [
            exact_or_glob(run_id, genome, "sqanti3/prefiltered/{run}_{genome}_prefiltered_classification.txt", "sqanti3/prefiltered/*_classification.txt"),
            exact_or_glob(run_id, genome, "sqanti3/postfiltered/{run}_{genome}_postfiltered_classification.txt", "sqanti3/postfiltered/*_classification.txt"),
        ])
    if stage == "dashboard":
        direct = status_for_patterns(genome_dir, ["dashboard/index.html", "dashboard/dashboard_data.json"])
        if direct == "complete":
            return direct
        return status_for_patterns(outdir, [f"dashboard/{genome}/index.html", f"dashboard/{genome}/dashboard_data.json"])
    raise ValueError(stage)


def recommended_common(statuses):
    recommendation = "none"
    for stage in STAGES:
        if stage == "bam_to_gtf":
            continue
        if all(stage_map.get(stage) == "complete" for stage_map in statuses.values()):
            recommendation = stage
    return recommendation


def recommended_for_genome(stage_map):
    recommendation = "none"
    for stage in STAGES:
        if stage == "bam_to_gtf":
            continue
        if stage_map.get(stage) == "complete":
            recommendation = stage
    return recommendation


def prepared_fasta_files(genome_dir):
    candidates = []
    for rel in [
        "reads/prepared_fasta",
        "reads/prepared",
        "reads/fasta",
        "reads/renamed",
        "reads/all_reads/fasta",
        "reads/all_reads/renamed",
    ]:
        candidates.extend(sorted((genome_dir / rel).glob("*.fasta.gz")))
    return [str(path) for path in candidates if path.is_file() and path.stat().st_size > 0]


def canonical_bams(genome_dir, run_id):
    bam_dir = genome_dir / "reads" / "canonical" / "chromosomes"
    bams = sorted(bam_dir.glob(f"{run_id}_*.canonical.bam"))
    chromosome_count = str(len(bams))
    output = []
    for bam in bams:
        bai = Path(str(bam) + ".csi")
        if not bai.exists():
            continue
        chromosome_id = bam.name.replace(f"{run_id}_", "").replace(".canonical.bam", "")
        output.append([str(bam), str(bai), chromosome_id, chromosome_count])
    return output


def chromosome_models(genome_dir, run_id):
    model_dir = genome_dir / "all_models" / "chromosomes"
    gtfs = sorted(model_dir.glob(f"{run_id}_*.combined.gtf"))
    chromosome_count = str(len(gtfs))
    output = []
    for gtf in gtfs:
        chromosome_id = gtf.name.replace(f"{run_id}_", "").replace(".combined.gtf", "")
        contains = model_dir / f"{run_id}_{chromosome_id}.combined.contains.tsv"
        if contains.exists():
            output.append([str(gtf), str(contains), chromosome_id, chromosome_count])
    return output


def full_length_chromosome_models(genome_dir, run_id):
    model_dir = genome_dir / "full_length" / "models" / "chromosomes"
    gtfs = sorted(model_dir.glob(f"{run_id}_*.filtered.gtf"))
    chromosome_count = str(len(gtfs))
    output = []
    for gtf in gtfs:
        chromosome_id = gtf.name.replace(f"{run_id}_", "")
        chromosome_id = chromosome_id.replace(".combined.filtered.gtf", "").replace(".filtered.gtf", "")
        contains = model_dir / f"{run_id}_{chromosome_id}_FL_contains.tsv"
        if contains.exists():
            output.append([str(gtf), str(contains), chromosome_id, chromosome_count])
    return output


def dashboard_bam_stats_inputs(genome_dir, run_id, genome_id):
    rows = []
    add_bam_stat_row(rows, "merged", genome_id, genome_id, genome_dir / "reads" / "all_reads" / "merged" / f"{run_id}_{genome_id}.bam")
    add_bam_stat_row(rows, "canonical_merged", genome_id, genome_id, genome_dir / "reads" / "canonical" / "merged" / f"{run_id}_{genome_id}.canonical.bam")
    add_bam_stat_row(rows, "full_length_merged", genome_id, genome_id, genome_dir / "full_length" / "reads" / "merged" / f"{run_id}_{genome_id}_full_length.bam")
    for bam in sorted((genome_dir / "reads" / "all_reads" / "split").glob("*.bam")):
        add_bam_stat_row(rows, "sample", genome_id, bam.stem, bam)
    for bam in sorted((genome_dir / "reads" / "all_reads" / "chromosomes").glob("*.bam")):
        add_bam_stat_row(rows, "chromosome", genome_id, bam.stem, bam)
    for bam in sorted((genome_dir / "reads" / "canonical" / "chromosomes").glob("*.canonical.bam")):
        chromosome_id = bam.name.replace(f"{run_id}_", "").replace(".canonical.bam", "")
        add_bam_stat_row(rows, "canonical", genome_id, chromosome_id, bam)
    for bam in sorted((genome_dir / "full_length" / "reads" / "chromosomes").glob("*.full_length.bam")):
        chromosome_id = bam.name.replace(f"{run_id}_", "").replace(".full_length.bam", "")
        add_bam_stat_row(rows, "full_length", genome_id, chromosome_id, bam)
    return rows


def add_bam_stat_row(rows, stage, genome_id, label, bam):
    if not bam.is_file() or bam.stat().st_size == 0:
        return
    bai = Path(str(bam) + ".csi")
    if not bai.is_file() or bai.stat().st_size == 0:
        raise SystemExit(f"Missing required BAM index for checkpoint dashboard stats: {bai}")
    rows.append([stage, genome_id, label, str(bam), str(bai)])


def exact_file(directory, filename):
    path = directory / filename
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"Missing required checkpoint file: {path}")
    return str(path)


def exact_or_glob(run_id, genome, exact_template, fallback):
    return exact_template.format(run=run_id, genome=genome) if run_id else fallback


def status_for_patterns(base, patterns):
    matches = [nonempty_matches(base, pattern) for pattern in patterns]
    if all(matches):
        return "complete"
    if any(path_hint_exists(base, pattern) or match for pattern, match in zip(patterns, matches)):
        return "partial"
    return "missing"


def nonempty_matches(base, pattern):
    return [path for path in base.glob(pattern) if path.is_file() and path.stat().st_size > 0]


def first_static_part(pattern):
    parts = []
    for part in Path(pattern).parts:
        if any(char in part for char in "*?["):
            break
        parts.append(part)
    return Path(*parts) if parts else Path(".")


def path_hint_exists(base, pattern):
    hint = base / first_static_part(pattern)
    return hint.exists() or (hint != base and hint.parent.exists())


def read_dicts(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def read_reference_sheet(path):
    with open(path, newline="") as handle:
        rows = [[value.strip() for value in row] for row in csv.reader(handle) if row]
    if not rows:
        return []
    first = [value.lower() for value in rows[0]]
    required = ["reference_file", "genome_id", "map_key", "species"]
    optional = ["reference_gtf", "polya_peak", "cage_peak"]
    if set(required).issubset(first):
        return [dict(zip(first, row)) for row in rows[1:]]
    columns_by_len = {
        4: required,
        5: required + ["reference_gtf"],
        6: required + ["reference_gtf", "polya_peak"],
        7: required + optional,
    }
    columns = columns_by_len.get(len(rows[0]))
    if not columns:
        raise SystemExit(f"Reference sheet {path} must have a recognised header or 4-7 headerless columns")
    return [dict(zip(columns, row)) for row in rows]


def group_metadata_by_genome(rows):
    grouped = {}
    for row in rows:
        grouped.setdefault(row.get("genome_id", ""), []).append(row)
    return grouped


def reference_lookup(references):
    output = {}
    for row in references:
        genome = row.get("genome_id", "")
        reference = row.get("reference_file", "")
        if genome and reference:
            output[genome] = reference
    return output


def seq_source_lookup(rows):
    lookup = {}
    for row in rows:
        reads = row.get("reads", "")
        if reads:
            lookup[Path(reads).name] = row.get("seq_source", "nanopore")
    return lookup


def seq_source_for_fasta(path, lookup, default):
    name = Path(path).name
    for sample_name, seq_source in lookup.items():
        stem = sample_name.replace(".fastq.gz", "").replace(".fq.gz", "").replace(".fasta.gz", "").replace(".fa.gz", "")
        if stem and stem in name:
            return seq_source
    return default


def first_value(rows, column, default):
    for row in rows:
        value = row.get(column, "")
        if value and value.lower() not in EMPTY_VALUES:
            return value
    return default


def write_rows(path, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerows(rows)


if __name__ == "__main__":
    main()
