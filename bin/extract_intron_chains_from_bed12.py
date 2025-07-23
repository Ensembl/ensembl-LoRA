#!/usr/bin/env python 
import sys
import pandas as pd
import argparse

def parse_bed12_line(line):
    fields = line.strip().split()
    chrom = fields[0]
    chrom_start = int(fields[1])
    chrom_end = int(fields[2])
    name = fields[3]
    score = fields[4]
    strand = fields[5]
    thick_start = int(fields[6])
    thick_end = int(fields[7])
    item_rgb = fields[8]
    block_count = int(fields[9])
    block_sizes = list(map(int, fields[10].split(',')))
    block_starts = list(map(int, fields[11].split(',')))
    
    return (chrom, chrom_start, chrom_end, name, score, strand, thick_start, thick_end, item_rgb, block_count, block_sizes, block_starts)

def extract_introns(chrom, chrom_start, block_count, block_sizes, block_starts):
    introns = []
    for i in range(block_count - 1):
        intron_start = chrom_start + block_starts[i] + block_sizes[i]
        intron_end = chrom_start + block_starts[i + 1]
        introns.append((intron_start, intron_end))
    return introns

def main():
    args = get_args()
    single_exon_data = []
    multi_exon_data = []
    with open(args.bed12, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            chrom, chrom_start, chrom_end, name, score, strand, thick_start, thick_end, item_rgb, block_count, block_sizes, block_starts = parse_bed12_line(line)
            if block_count > 1:
                introns = extract_introns(chrom, chrom_start, block_count, block_sizes, block_starts)
                if strand == '-':
                    introns.reverse()
                intron_chain = ','.join([f'{intron_start + 1}-{intron_end}' for intron_start, intron_end in introns])
                multi_exon_data.append([chrom, name, f'({strand})({chrom}){intron_chain}'])
            else:
                exon_start = chrom_start + block_starts[0]
                exon_end = exon_start + block_sizes[0]
                if strand == '+':
                    exon_start, exon_end = exon_start + 1, exon_end
                    single_exon_data.append([chrom, name, f'({strand})({chrom}){exon_start}-{exon_end}'])
                if strand == '-':
                    exon_start, exon_end = exon_end, exon_start + 1
                    single_exon_data.append([chrom, name, f'({strand})({chrom}){exon_start}-{exon_end}'])

    single_exon_df = pd.DataFrame(single_exon_data, columns=['seqname', 'transcript_id', 'intron_chain'])
    multi_exon_df = pd.DataFrame(multi_exon_data, columns=['seqname', 'transcript_id', 'intron_chain'])
    single_exon_df.to_csv(f'{args.prefix}_single_exons.tsv', sep="\t", index=False, header=False)
    multi_exon_df.to_csv(f'{args.prefix}_intron_chains.tsv', sep="\t", index=False, header=False)

def get_args():
    parser = argparse.ArgumentParser(description="This script takes as input the reads sample sheet and the reference sample sheet and outputs the combined sample sheet")
    parser.add_argument("bed12", help="The bed12 for the reads")
    parser.add_argument("prefix", help="Output Prefix")
    return parser.parse_args()

if __name__ == '__main__':
    main()