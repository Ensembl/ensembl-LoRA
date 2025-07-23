#!/usr/bin/env python
import pandas as pd
import os
import argparse

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process picard output')
    parser.add_argument('genome_id', type=str, help='input file')
    parser.add_argument('individual_bin_count', type=str, help='output file')
    parser.add_argument('total_bin_count', type=str, help='output file')
    parser.add_argument('output', type=str, help='output file')
    args = parser.parse_args()
    
    # Assign command line arguments to variables
    genome_id = args.genome_id
    split_files = [os.path.join(os.getcwd(), f) for f in os.listdir() if f.endswith('.fixed.gtf')]
    bin_count = args.total_bin_count
    output = args.output
    
    # Create DataFrame
    data = []
    for split_file in split_files:
        data.append({
            'binned_gtf': split_file,
            'genome_id': genome_id,
            'individual_bin_count': args.individual_bin_count,
            'total_bin_count': bin_count,
        })
    picard_df = pd.DataFrame(data)
    picard_df.to_csv(output, index=False, header=False)

if __name__ == '__main__':
    main()