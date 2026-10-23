#!/usr/bin/env python
import argparse
import subprocess
import sys


def main():
    args = get_args()
    strand_info = load_strand_info()
    write_canonical_bam(args.bam, args.prefix, strand_info, args.threads)
    run(["samtools", "index", "-c", "-@", str(args.threads), f"{args.prefix}.canonical.bam"])


def get_args():
    parser = argparse.ArgumentParser(description="Fix BAM strand flags for reads with canonical splice-site evidence")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--threads", type=int, default=1)
    return parser.parse_args()


def load_strand_info():
    strand_info = {}
    fixed_beds = sorted(path for path in subprocess.check_output(["find", ".", "-maxdepth", "1", "-name", "*.fixed.bed"], text=True).splitlines() if path)
    if not fixed_beds:
        return strand_info
    for bed_path in fixed_beds:
        with open(bed_path, encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    continue
                strand_info[fields[3]] = fields[5]
    return strand_info


def write_canonical_bam(bam, prefix, strand_info, threads):
    header = subprocess.run(["samtools", "view", "-H", bam], check=True, stdout=subprocess.PIPE, text=True).stdout
    samtools_out = subprocess.Popen(
        ["samtools", "view", "-@", str(threads), "-b", "-o", f"{prefix}.canonical.bam", "-"],
        stdin=subprocess.PIPE,
        text=True,
    )
    assert samtools_out.stdin is not None
    with samtools_out.stdin as out_handle:
        for line in header.splitlines():
            if not line.startswith("@PG"):
                out_handle.write(line + "\n")
        samtools_in = subprocess.Popen(["samtools", "view", bam], stdout=subprocess.PIPE, text=True)
        assert samtools_in.stdout is not None
        for line in samtools_in.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            read_id = fields[0]
            strand = strand_info.get(read_id)
            if strand is None:
                continue
            flag = int(fields[1])
            if strand == "+":
                flag &= ~16
            elif strand == "-":
                flag |= 16
            fields[1] = str(flag)
            out_handle.write("\t".join(fields) + "\n")
        samtools_in.stdout.close()
        if samtools_in.wait() != 0:
            raise SystemExit(f"samtools view failed while reading {bam}")
    if samtools_out.wait() != 0:
        raise SystemExit(f"samtools view failed while writing {prefix}.canonical.bam")


def run(command):
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"ERROR: command failed: {' '.join(command)}", file=sys.stderr)
        raise SystemExit(exc.returncode)


if __name__ == "__main__":
    main()
