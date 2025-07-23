#!/usr/bin/env python 
import pandas as pd
import argparse
import os

def main():
    options = getOptions()
    
    try:
        # Extracting chromosome ids from GTF file
        chromosome_ids = extract_chromosome_ids(options.gtf_file)
        
        # Creating a dataframe
        chromosome_df = pd.DataFrame(chromosome_ids)
        chromosome_df['genome_id'] = options.genome_id
        if os.path.exists(options.gtf_file):
            chromosome_df['gtf'] = os.path.join(os.getcwd(), options.gtf_file)
        else:
            raise FileNotFoundError(f"GTF file not found: {options.gtf_file}")
        if os.path.exists(options.contains_file):
            chromosome_df['contains'] = os.path.join(os.getcwd(), options.contains_file)
        else:
            raise FileNotFoundError(f"Contains file not found: {options.contains_file}")
        chromosome_df['chromosome_count'] = len(chromosome_ids)
        chromosome_df['chromosome_count'] = chromosome_df['chromosome_count'].astype(int)
        chromosome_df.columns = ['chromosome', 'genome_id', 'gtf', 'contains', 'chromosome_count']
        
        # Saving the dataframe to a CSV file
        chromosome_df.to_csv(options.output, index=False, header = False)
    except FileNotFoundError as fnf_error:
        print(f"File not found: {fnf_error}")
    except pd.errors.EmptyDataError as ede_error:
        print(f"Empty data error: {ede_error}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return

def getOptions():
    parser = argparse.ArgumentParser(description="Extract chromosome ids from a GTF file and create CSV file for nextflow")
    parser.add_argument("--genome_id", "-g", type=str, help="Genome ID", required=True)
    parser.add_argument("--gtf_file", "-f", type=str, help="GTF file", required=True)
    parser.add_argument("--contains_file", "-c", type=str, help="Contains file", required=True)
    parser.add_argument("--output", "-o", type=str, help="Output file", required=True)
    
    options = parser.parse_args()
    return options

def extract_chromosome_ids(gtf_file):
    chromosome_ids = set()
    with open(gtf_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                columns = line.split('\t')
                chromosome_ids.add(columns[0])
    return list(chromosome_ids)

if __name__ == "__main__":
    main()