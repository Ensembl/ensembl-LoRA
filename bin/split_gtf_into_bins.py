#!/usr/bin/env python
import os
import argparse
from collections import defaultdict

def split_gtf_into_bins(gtf_file, num_bins, run_id):
    transcripts = defaultdict(list)
    # Read the GTF file and group exons by transcript
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            attributes = fields[8]
            transcript_id = None
            for attr in attributes.split(';'):
                if 'transcript_id' in attr:
                    transcript_id = attr.split('"')[1]
                    break
            if transcript_id:
                transcripts[transcript_id].append(line)
    # Calculate the number of transcripts per bin
    total_transcripts = len(transcripts)
    transcripts_per_bin = total_transcripts // num_bins
    if total_transcripts % num_bins != 0:
        transcripts_per_bin += 1
    # Split transcripts into bins
    bin_count = 0
    bin_lines = []
    current_bin_size = 0
    for transcript_id, lines in transcripts.items():
        if current_bin_size >= transcripts_per_bin:
            bin_count += 1
            output_file = f'{run_id}_{bin_count:04d}.fixed.gtf'
            with open(output_file, 'w') as out_f:
                out_f.writelines(bin_lines)
            bin_lines = []
            current_bin_size = 0
        bin_lines.extend(lines)
        current_bin_size += 1
    # Write remaining lines to the last bin if any
    if bin_lines:
        bin_count += 1
        output_file = f'{run_id}_{bin_count:04d}.fixed.gtf'
        with open(output_file, 'w') as out_f:
            out_f.writelines(bin_lines)

def main():
    parser = argparse.ArgumentParser(description='Split GTF file into bins')
    parser.add_argument('gtf_file', type=str, help='GTF file to split')
    parser.add_argument('num_bins', type=int, help='Number of bins')
    parser.add_argument('run_id', type=str, help='Run ID')
    args = parser.parse_args()
    split_gtf_into_bins(args.gtf_file, args.num_bins, args.run_id)

if __name__ == "__main__":
    main()