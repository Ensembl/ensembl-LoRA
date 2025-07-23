#!/usr/bin/env python 
import pandas as pd
import numpy as np
import argparse
import os

def main():
    # Get arguments
    args = get_args()
    
    # Get chromosomes and read counts using samtools idxstats
    os.system(f'samtools idxstats {args.bam} > {args.prefix}_idxstats.txt')
    idx_stats = pd.read_csv(f'{args.prefix}_idxstats.txt', sep="\t", header=None, names=["chromosome", "length", "mapped_reads", "unmapped_reads"])
    
    # Filter out empty chromosomes
    idx_stats = idx_stats[idx_stats['mapped_reads'] > 0]
    
    # Divide reads by bin size
    idx_stats['bin_count'] = (idx_stats['mapped_reads'] / int(args.bin_size)).apply(np.ceil).astype(int)
    
    # Add information of idxstats
    idx_stats['bam'] = os.path.join(os.getcwd(),args.bam)
    idx_stats['bai'] = os.path.join(os.getcwd(),args.bai)
    idx_stats['reference'] = os.path.join(os.getcwd(),args.reference_file)
    idx_stats['genome_id'] = args.genome_id
    idx_stats['chromosome_count'] = idx_stats.shape[0]
    idx_stats['chromosome_group'] = idx_stats['genome_id'].astype(str) + '_' + idx_stats['chromosome'].astype(str)
    
    # Reorder and select relevant columns
    idx_stats = idx_stats[['bam', 'bai', 'reference', 'genome_id', 'chromosome', 'chromosome_group', 'chromosome_count', 'bin_count']]
    
    # Save sample sheet
    idx_stats.to_csv(f'{args.prefix}_sample_sheet.csv', index=False, header=False)

def get_args():
    parser = argparse.ArgumentParser(description="This script creates a sample sheet containing chromosome and bin count information for a merged bam file")
    parser.add_argument("genome_id", help="Genome ID")
    parser.add_argument("bam", help="Bam file")
    parser.add_argument("bai", help="Bam index file")
    parser.add_argument("reference_file", help="Reference file")
    parser.add_argument("bin_size", help="Bin size")
    parser.add_argument("prefix", help="Output prefix")
    return parser.parse_args()

if __name__ == "__main__":
    main()

