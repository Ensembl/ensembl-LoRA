#!/usr/bin/env python 
import pandas as pd
import argparse

def main():
    # Processing the arguments
    args = get_args()
    
    # Reading the sample sheets
    reads_sample_sheet = pd.read_csv(args.reads_sample_sheet, header=None, names=['reads', 'format', 'seq_source', 'map_key', 'species'])
    reference_sample_sheet = pd.read_csv(args.reference_sample_sheet, header=None, names=['reference_file','genome_id', 'map_key', 'species'])
    
    # Running checks on the reads sample sheet
    reads_sample_sheet_checks(reads_sample_sheet, reference_sample_sheet)
    reference_sample_sheet_checks(reference_sample_sheet, reads_sample_sheet)
    
    # Combining the sample sheets
    merged_sample_sheets = pd.merge(reads_sample_sheet, reference_sample_sheet, on=['map_key', 'species'], how='inner')
    
    # Counting number of files per genome_id
    genome_id_counts = merged_sample_sheets.groupby('genome_id').size().reset_index(name='counts')
    merged_sample_sheets = pd.merge(merged_sample_sheets, genome_id_counts, on='genome_id', how='inner')
    
    # Ordering required columns
    merged_sample_sheets = merged_sample_sheets[['reads', 'format', 'seq_source', 'reference_file', 'genome_id', 'counts', 'species']]
    
    # Writing the output sample sheet
    merged_sample_sheets.to_csv(args.output_sample_sheet, header=False, index=False, sep=',')

def get_args():
    parser = argparse.ArgumentParser(description="This script takes as input the reads sample sheet and the reference sample sheet and outputs the combined sample sheet")
    parser.add_argument("reads_sample_sheet", help="The sample sheet for the reads")
    parser.add_argument("reference_sample_sheet", help="The sample sheet for the reference")
    parser.add_argument("output_sample_sheet", help="The output sample sheet")
    return parser.parse_args()

def reads_sample_sheet_checks(sample_sheet, reference_sample_sheet):
    # Checking if the sample sheet has the required columns
    if 'reads' not in sample_sheet.columns or 'format' not in sample_sheet.columns or 'seq_source' not in sample_sheet.columns or 'map_key' not in sample_sheet.columns or 'species' not in sample_sheet.columns:
        raise ValueError("The reads sample sheet is missing one of the required columns: reads, format, seq_source, map_key, species")
    
    # Check that the reads files exist and are not corrupt
    for reads_file in sample_sheet['reads']:
        try:
            with open(reads_file) as f:
                pass
        except FileNotFoundError:
            raise FileNotFoundError(f"The reads file {reads_file} does not exist")
        except Exception as e:
            raise Exception(f"The reads file {reads_file} is corrupt")
    
    # Checking if the required columns have any missing values
    if sample_sheet['reads'].isnull().values.any() or sample_sheet['format'].isnull().values.any() or sample_sheet['seq_source'].isnull().values.any() or sample_sheet['map_key'].isnull().values.any() or sample_sheet['species'].isnull().values.any():
        raise ValueError("The reads sample sheet has missing values in one of the required columns: reads, format, seq_source, map_key, species")
    
    # Check that the seq_source column is either 'pacbio' or 'nanopore'
    if not sample_sheet['seq_source'].isin(['pacbio', 'nanopore']).all():
        raise ValueError("The seq_source column in the reads sample sheet should have values either 'pacbio' or 'nanopore'")
    
    # Check that the provided genome_id exists in the reference sample sheet
    if not sample_sheet['map_key'].isin(reference_sample_sheet['map_key']).all():
        raise ValueError("The genome_id in the reads sample sheet should exist in the reference sample sheet")
    
    # Check that the provided species exists in the reference sample sheet
    if not sample_sheet['species'].isin(reference_sample_sheet['species']).all():
        raise ValueError("The map_key in the reads sample sheet should exist in the reference sample sheet")

def reference_sample_sheet_checks(reference_sample_sheet, sample_sheet):
    # Checking if the reference sample sheet has the required columns
    if 'reference_file' not in reference_sample_sheet.columns or 'genome_id' not in reference_sample_sheet.columns or 'map_key' not in reference_sample_sheet.columns or 'species' not in reference_sample_sheet.columns:
        raise ValueError("The reference sample sheet is missing one of the required columns: reference_file, genome_id, map_key, species")

    # Check that the reference files exist and are not corrupt
    for reference_file in reference_sample_sheet['reference_file']:
        try:
            with open(reference_file) as f:
                pass
        except FileNotFoundError:
            raise FileNotFoundError(f"The reference file {reference_file} does not exist")
        except Exception as e:
            raise Exception(f"The reference file {reference_file} is corrupt")

    # Checking if the required columns have any missing values
    if reference_sample_sheet['reference_file'].isnull().values.any() or reference_sample_sheet['genome_id'].isnull().values.any() or reference_sample_sheet['map_key'].isnull().values.any() or reference_sample_sheet['species'].isnull().values.any():
        raise ValueError("The reference sample sheet has missing values in one of the required columns: reference_file, genome_id, map_key, species")

if __name__ == "__main__":
    main()