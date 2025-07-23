#!/usr/bin/env python
import pandas as pd
import argparse

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Merge TMERGE contains outputs")
    parser.add_argument("tmerge", help="Initial TMERGE contains file")
    parser.add_argument("retmerge", help="Final TMERGE contains file")
    parser.add_argument("prefix", help="Prefix for output")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Assing the aguments to variables
    initial = args.tmerge
    final = args.retmerge
    prefix = args.prefix
    
    # Importing contains columns
    initial_contains = pd.read_csv(initial, sep = "\t", header = None)
    initial_contains.columns = ['initial_chromosome', 'initial_transcript_id', 'reads']
    
    final_contains = pd.read_csv(final, sep = "\t", header = None)
    final_contains.columns = ['chromosome', 'transcript_id', 'contains']
    
    # Pivot rows with commas
    rows_with_comma = final_contains[final_contains['contains'].str.contains(',')]
    rows_without_comma = final_contains[~final_contains['contains'].str.contains(',')]
    rows_with_comma = final_contains.assign(contains=final_contains['contains'].str.split(',')).explode('contains')
    
    # Combine the two dataframes
    final_contains_long = pd.concat([rows_with_comma, rows_without_comma])
    
    # Merging the two tables 
    merged_table = pd.merge(final_contains_long, initial_contains, left_on=["chromosome", 'contains'], right_on=["initial_chromosome",'initial_transcript_id'], how = 'inner')
    merged_table = merged_table[['chromosome', 'transcript_id', 'reads']].drop_duplicates()
    
    # Combining reads for each transcript
    grouped_df = merged_table.groupby(["chromosome", 'transcript_id'])['reads'].agg(lambda x: ','.join(map(str, x))).reset_index()
    
    # Saving to file
    grouped_df.to_csv(f'{prefix}.contains.tsv', sep = "\t", index = False, header=False)
    
if __name__ == "__main__":
    main()