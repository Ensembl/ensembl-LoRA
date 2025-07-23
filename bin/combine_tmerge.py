#!/usr/bin/env python
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Process some GTF and TSV files.')
    # parser.add_argument('--gtf', '-g', type=str, nargs='+', help='Path to the GTF file')
    # parser.add_argument('--contains', '-c', type=str, nargs='+', help='Path to the contains TSV file')
    parser.add_argument('--run_name', '-r', type=str, help='Run name (first or other)')
    parser.add_argument('--output_prefix', '-o', type=str, help='Output prefix for the files')

    args = parser.parse_args()

    if args.run_name == "first":
        gtf_name = f"{args.output_prefix}_initial_transcripts.gtf"
        contains_name = f"{args.output_prefix}_initial_contains.tsv"
    elif args.run_name == "filter":
        gtf_name = f"{args.output_prefix}_filtered_transcripts.gtf"
        contains_name = f"{args.output_prefix}_filtered_contains.tsv"
    elif args.run_name == "final":
        gtf_name = f"{args.output_prefix}.gtf"
        contains_name = f"{args.output_prefix}_contains.tsv"
    else:
        gtf_name = f"{args.output_prefix}_unfiltered.gtf"
        contains_name = f"{args.output_prefix}_unfiltered_contains.tsv"
    if args.run_name:
        os.system(f"cat *gtf | sort -k1,1V -k4,4n > {gtf_name}")
    os.system(f"cat *tsv | awk -F'\\t' '{{split($2, a, \"_\");split(a[length(a)], b, \".\");print b[1] \"\\t\" b[2] \"\\t\" $0}}' | sort -k3,3V -k1,1n -k2,2n | cut -f3- > {contains_name}")

if __name__ == "__main__":
    main()
