#!/usr/bin/env python 
import pandas as pd
import re
import argparse

def main():
    args = get_args()
    # Read in the files
    models_intron_chains_contain = pd.read_csv(args.models_intron_chains, sep='\t', header=None, names=['chromosome', 'model','intron_chain'], dtype=str)
    models_single_exons_contain = pd.read_csv(args.models_single_exons, sep='\t', header=None, names=['chromosome', 'model','intron_chain'], dtype=str)
    reads_intron_chains_contain = pd.read_csv(args.reads_intron_chains, sep='\t', header=None, names=['chromosome', 'read','intron_chain'], dtype=str)
    reads_single_exons_contain = pd.read_csv(args.reads_single_exons, sep='\t', header=None, names=['chromosome', 'read','intron_chain'], dtype=str)
    all_contains = pd.read_csv(args.all_contains, sep='\t', header=None, names=['chromosome', 'model','contains'], dtype=str)
    
    # Pivot rows with commas
    rows_with_comma = all_contains[all_contains['contains'].str.contains(',')]
    rows_without_comma = all_contains[~all_contains['contains'].str.contains(',')]
    rows_with_comma = rows_with_comma.assign(contains=rows_with_comma['contains'].str.split(',')).explode('contains')
    
    # Combine the two dataframes
    final_contains_long = pd.concat([rows_with_comma, rows_without_comma])
    final_contains_long.columns = ['chromosome', 'model', 'read']
    
    # Add intron chains to reads
    final_contains_long = final_contains_long.merge(reads_intron_chains_contain, on=['chromosome', 'read'], how='left')
    
    # Split into single exon and intron chain
    final_contains_intron_chain = final_contains_long[~final_contains_long['intron_chain'].isnull()]
    final_contains_single_exon = final_contains_long[final_contains_long['intron_chain'].isnull()]
    final_contains_single_exon = final_contains_single_exon.drop(columns=['intron_chain'])
    final_contains_single_exon = final_contains_single_exon.merge(reads_single_exons_contain, on=['chromosome', 'read'], how='left')
    
    # Initial stats
    total_reads = final_contains_long['read'].nunique()
    total_multi_exon = final_contains_intron_chain['read'].nunique()
    total_single_exon = final_contains_single_exon['read'].nunique()
    
    # Extract full length reads
    final_contains_intron_chain = final_contains_intron_chain.merge(models_intron_chains_contain, on=['chromosome', 'model', 'intron_chain'], how='inner')
    final_contains_single_exon = final_contains_single_exon.merge(models_single_exons_contain, on=['chromosome', 'model', 'intron_chain'], how='inner')
    
    # Second stats
    total_full_length_multi_exon = final_contains_intron_chain['read'].nunique()
    total_full_length_single_exon = final_contains_single_exon['read'].nunique()
    total_full_length = total_full_length_multi_exon + total_full_length_single_exon
    
    # Extract transcript ids with only 1 full length read
    full_length_multi_exon_read_counts = final_contains_intron_chain.groupby('model')['read'].nunique().reset_index()
    full_length_single_exon_read_counts = final_contains_single_exon.groupby('model')['read'].nunique().reset_index()
    
    # Third stats
    full_length_multi_exon_multi_read = full_length_multi_exon_read_counts[full_length_multi_exon_read_counts['read'] > 1]
    full_length_single_exon_multi_read = full_length_single_exon_read_counts[full_length_single_exon_read_counts['read'] > 1]
    
    # Output stats
    pd.DataFrame({
        'Sample': [args.output_prefix],
        'Total Reads': [total_reads],
        'Total Multi-Exon Reads': [total_multi_exon],
        'Total Single-Exon Reads': [total_single_exon],
        'Total Full Length Reads': [total_full_length],
        'Total Full Length Multi-Exon Reads': [total_full_length_multi_exon],
        'Total Full Length Single-Exon Reads': [total_full_length_single_exon],
        'Total Transcripts': [all_contains['model'].nunique()],
        'Total Multi-Exon Transcripts with more that 1 FL read': [full_length_multi_exon_multi_read.shape[0]],
        'Total Single-Exon Transcripts with more that 1 FL read': [full_length_single_exon_multi_read.shape[0]]
    }).to_csv(f'{args.output_prefix}_stats.tsv', sep='\t', index=False)

    # Create full length reads contains
    full_length_contains = pd.concat([final_contains_intron_chain, final_contains_single_exon]).reset_index(drop=True)
    full_length_contains = full_length_contains.drop(columns=['intron_chain'])
    full_length_contains = full_length_contains.groupby(["chromosome", 'model'])['read'].agg(lambda x: ','.join(map(str, x))).reset_index()

    # Sort by chromosome and model number
    full_length_contains['model_number'] = full_length_contains['model'].apply(lambda x: int(re.search(r'_(\d+)$', x).group(1)))
    full_length_contains = full_length_contains.sort_values(['chromosome', 'model_number'])
    full_length_contains.drop(columns=['model_number'], inplace=True)

    # Exporting full length reads contains
    full_length_contains.to_csv(f'{args.output_prefix}_full_length_contains.tsv', sep='\t', index=False, header=False)

def get_args():
    parser = argparse.ArgumentParser(description="This script takes in the output of the merging transcripts and outputs the stats and full length reads contains")
    parser.add_argument("models_intron_chains", help="Intron chains for models")
    parser.add_argument("models_single_exons", help="Single exon models")
    parser.add_argument("reads_intron_chains", help="Intron chains for reads")
    parser.add_argument("reads_single_exons", help="Single exon reads")
    parser.add_argument("all_contains", help="Contains file")
    parser.add_argument("output_prefix", help="Output prefix")
    
    return parser.parse_args()

if __name__ == '__main__':
    main()

