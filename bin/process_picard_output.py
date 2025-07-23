#!/usr/bin/env python
import pandas as pd
import os
import argparse

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process picard output')
    parser.add_argument('genome_id', type=str, help='input file')
    parser.add_argument('reference_genome', type=str, help='reference genome')
    parser.add_argument('chromosome_id', type=str, help='crhomosome id')
    parser.add_argument('chromosome_count', type=str, help='chromosome count')
    parser.add_argument('bin_count', type=int, help='bin count')
    parser.add_argument('output', type=str, help='output file')
    args = parser.parse_args()
    
    # Assign command line arguments to variables
    genome_id = args.genome_id
    picard_bams = [os.path.join(os.getcwd(), f) for f in os.listdir() if f.endswith('picard.bam')]
    reference_genome = os.path.abspath(args.reference_genome)
    chromosome_id = args.chromosome_id
    chromosome_count = args.chromosome_count
    bin_count = args.bin_count
    output = args.output
    
    # Create DataFrame
    data = []
    for bam in picard_bams:
        data.append({
            'picard_bam': bam,
            'reference_genome': reference_genome,
            'genome_id': genome_id,
            'chromosome_id': chromosome_id,
            'chromosome_count': chromosome_count,
            'bin_count': bin_count
        })
    picard_df = pd.DataFrame(data)
    picard_df.to_csv(output, index=False, header=False)

if __name__ == '__main__':
    main()