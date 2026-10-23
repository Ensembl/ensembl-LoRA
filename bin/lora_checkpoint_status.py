#!/usr/bin/env python3
import argparse
import csv
import os
from pathlib import Path


STAGES = [
    "read_prep",
    "mapped_bam",
    "bam_to_gtf",
    "merged_models",
    "full_length",
    "sqanti3",
    "dashboard",
]
RESUMABLE_STAGES = {"read_prep", "mapped_bam", "merged_models", "full_length", "sqanti3", "dashboard"}

NON_GENOME_DIRS = {"dashboard", "pipeline_info", "summary", "work", "tmp", "checkpoints"}


def main():
    args = get_args()
    outdir = resolve_outdir(args.outdir)
    genomes = [args.genome] if args.genome else discover_genomes(outdir)
    if not genomes:
        raise SystemExit(no_genomes_message(outdir, args.run_id))

    statuses = {genome: {stage: stage_status(outdir, args.run_id, genome, stage) for stage in STAGES} for genome in genomes}
    common = recommended_common(statuses)
    per_genome = {genome: recommended_for_genome(stage_map) for genome, stage_map in statuses.items()}

    if args.format == "json":
        import json

        print(json.dumps({
            "outdir": str(outdir),
            "recommended_common_resume_from": common,
            "recommended_resume_from_per_genome": per_genome,
            "genomes": statuses,
        }, indent=2))
        return

    print(f"outdir\t{outdir}")
    print(f"recommended_common_resume_from\t{common}")
    print("recommended_resume_from_per_genome\t" + ",".join(f"{genome}:{stage}" for genome, stage in per_genome.items()))
    print()
    writer = csv.writer(os.sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(["genome", *STAGES, "recommended_resume_from"])
    for genome in genomes:
        writer.writerow([genome, *[statuses[genome][stage] for stage in STAGES], per_genome[genome]])


def get_args():
    parser = argparse.ArgumentParser(description="Audit LoRA published-output checkpoint completeness")
    parser.add_argument("--outdir", required=True, help="LoRA output directory")
    parser.add_argument("--run_id", default="", help="Optional run ID for exact filename checks")
    parser.add_argument("--genome", default="", help="Restrict audit to one genome")
    parser.add_argument("--format", choices=["table", "json"], default="table")
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
        path.name
        for path in outdir.iterdir()
        if path.is_dir() and path.name not in NON_GENOME_DIRS and not path.name.startswith(".")
    )


def stage_status(outdir, run_id, genome, stage):
    genome_dir = outdir / genome
    if stage == "read_prep":
        read_summary = status_for_patterns(genome_dir, [
            "reads/read_summaries/*_read_summary.tsv",
        ])
        prepared_fasta = "complete" if prepared_fasta_files(genome_dir) else "missing"
        if read_summary == "complete" and prepared_fasta == "complete":
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
        direct = status_for_patterns(genome_dir, [
            "dashboard/index.html",
            "dashboard/dashboard_data.json",
        ])
        if direct == "complete":
            return direct
        return status_for_patterns(outdir, [
            f"dashboard/{genome}/index.html",
            f"dashboard/{genome}/dashboard_data.json",
        ])
    raise ValueError(stage)


def exact_or_glob(run_id, genome, exact_template, fallback):
    if run_id:
        return exact_template.format(run=run_id, genome=genome)
    return fallback


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
    if hint.exists():
        return True
    if hint != base and hint.parent.exists():
        return True
    return False


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
    return [path for path in candidates if path.is_file() and path.stat().st_size > 0]


def recommended_common(statuses):
    recommendation = "none"
    for stage in STAGES:
        if stage not in RESUMABLE_STAGES:
            continue
        if all(stage_map.get(stage) == "complete" for stage_map in statuses.values()):
            recommendation = stage
    return recommendation


def recommended_for_genome(stage_map):
    recommendation = "none"
    for stage in STAGES:
        if stage not in RESUMABLE_STAGES:
            continue
        if stage_map.get(stage) == "complete":
            recommendation = stage
    return recommendation


if __name__ == "__main__":
    main()
