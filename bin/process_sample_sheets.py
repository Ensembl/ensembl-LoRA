#!/usr/bin/env python 
import pandas as pd
import argparse
import gzip
import csv
import os
import re

def main():
    # Processing the arguments
    args = get_args()
    
    # Reading the sample sheets
    reads_sample_sheet = read_reads_sample_sheet(args.reads_sample_sheet)
    reference_sample_sheet = read_reference_sample_sheet(args.reference_sample_sheet)
    
    # Running checks on the reads sample sheet
    reads_sample_sheet_checks(reads_sample_sheet, reference_sample_sheet)
    reference_sample_sheet_checks(reference_sample_sheet, reads_sample_sheet)
    
    # Combining the sample sheets
    merged_sample_sheets = pd.merge(reads_sample_sheet, reference_sample_sheet, on=['map_key', 'species'], how='inner')
    validate_sample_sheet_genome_ids(merged_sample_sheets)
    
    # Counting number of files per genome_id
    genome_id_counts = merged_sample_sheets.groupby('genome_id').size().reset_index(name='counts')
    merged_sample_sheets = pd.merge(merged_sample_sheets, genome_id_counts, on='genome_id', how='inner')
    
    # Ordering required columns. The final read_id_map field is optional and is
    # consumed only by FASTA/FASTQ read renaming.
    merged_sample_sheets = merged_sample_sheets[['reads', 'format', 'seq_source', 'reference_file', 'genome_id', 'counts', 'species', 'tissue', 'cell_line', 'read_id_map']]
    
    # Writing the output sample sheet
    merged_sample_sheets.to_csv(args.output_sample_sheet, header=False, index=False, sep=',')

    sample_metadata = pd.merge(reads_sample_sheet, reference_sample_sheet[['genome_id', 'map_key', 'species']], on=['map_key', 'species'], how='inner')
    validate_sample_sheet_genome_ids(sample_metadata)
    sample_metadata = sample_metadata[['reads', 'format', 'seq_source', 'map_key', 'species', 'genome_id', 'tissue', 'cell_line', 'sample_id', 'input_file_id', 'read_id_map']]
    sample_metadata.to_csv(args.sample_metadata, header=True, index=False, sep=',')

    sqanti_reference_sheet = reference_sample_sheet[
        reference_sample_sheet['reference_gtf'].apply(is_provided)
    ][['reference_file', 'genome_id', 'reference_gtf', 'polya_peak', 'cage_peak']].drop_duplicates(subset=['genome_id'])
    sqanti_reference_sheet.to_csv(args.sqanti_reference_sheet, header=False, index=False, sep=',')

def get_args():
    parser = argparse.ArgumentParser(description="This script takes as input the reads sample sheet and the reference sample sheet and outputs the combined sample sheet")
    parser.add_argument("reads_sample_sheet", help="The sample sheet for the reads")
    parser.add_argument("reference_sample_sheet", help="The sample sheet for the reference")
    parser.add_argument("output_sample_sheet", help="The output sample sheet")
    parser.add_argument("--sqanti_reference_sheet", default="sqanti3_reference_sheets.csv", help="The output reference sheet for SQANTI3-enabled references")
    parser.add_argument("--sample_metadata", default="sample_metadata.csv", help="The output sample metadata table for reporting")
    return parser.parse_args()

def read_reads_sample_sheet(reads_sample_sheet):
    required = ['reads', 'format', 'seq_source', 'map_key', 'species']
    optional = ['genome_id', 'tissue', 'cell_line', 'sample_id', 'read_id_map']
    sample_sheet = read_csv_with_optional_header(reads_sample_sheet, required, optional, {
        5: required,
        7: required + ['tissue', 'cell_line'],
        8: required + ['tissue', 'cell_line', 'read_id_map'],
        9: required + ['tissue', 'cell_line', 'sample_id', 'read_id_map'],
        10: required + optional,
    }, "reads sample sheet")
    sample_sheet = clean_string_columns(sample_sheet)
    sample_sheet['format'] = sample_sheet['format'].str.lower()
    sample_sheet['seq_source'] = sample_sheet['seq_source'].str.lower()
    for column in optional:
        if column not in sample_sheet.columns:
            sample_sheet[column] = 'NA'
        sample_sheet[column] = sample_sheet[column].fillna('NA')
    sample_sheet = add_sample_identity(sample_sheet)
    sample_sheet = sample_sheet.rename(columns={'genome_id': 'sample_sheet_genome_id'})
    return sample_sheet[required + ['sample_sheet_genome_id', 'tissue', 'cell_line', 'sample_id', 'input_file_id', 'read_id_map']]

def read_reference_sample_sheet(reference_sample_sheet):
    required = ['reference_file', 'genome_id', 'map_key', 'species']
    optional = ['reference_gtf', 'polya_peak', 'cage_peak']
    reference_sample_sheet = read_csv_with_optional_header(reference_sample_sheet, required, optional, {
        4: required,
        5: required + ['reference_gtf'],
        6: required + ['reference_gtf', 'polya_peak'],
        7: required + optional,
    }, "reference sample sheet")
    reference_sample_sheet = clean_string_columns(reference_sample_sheet)
    for column in optional:
        if column not in reference_sample_sheet.columns:
            reference_sample_sheet[column] = 'NA'
        reference_sample_sheet[column] = reference_sample_sheet[column].fillna('NA')
    reference_sample_sheet = reference_sample_sheet[required + optional]
    return reference_sample_sheet

def read_csv_with_optional_header(path, required_columns, optional_columns, headerless_column_sets, sheet_name):
    with open(path, newline='') as handle:
        rows = [[str(value).strip() for value in row] for row in csv.reader(handle)]
    rows = [row for row in rows if row and any(value != '' for value in row)]
    if not rows:
        raise ValueError(f"The {sheet_name} is empty")

    first_row = [normalize_header_name(value) for value in rows[0]]
    expected_columns = set(required_columns + optional_columns)
    has_header = set(required_columns).issubset(first_row) and set(first_row).issubset(expected_columns)

    if has_header:
        header = first_row
        missing = [column for column in required_columns if column not in header]
        if missing:
            raise ValueError(f"The {sheet_name} is missing required columns: {', '.join(missing)}")
        unknown = [column for column in header if column not in expected_columns]
        if unknown:
            raise ValueError(f"The {sheet_name} has unknown columns: {', '.join(unknown)}")
        bad_rows = [str(index) for index, row in enumerate(rows[1:], start=2) if len(row) != len(header)]
        if bad_rows:
            raise ValueError(f"The {sheet_name} has rows with a different number of columns from the header: {', '.join(bad_rows)}")
        sheet = pd.DataFrame(rows[1:], columns=header)
        return sheet

    row_lengths = {len(row) for row in rows}
    if len(row_lengths) != 1:
        raise ValueError(f"The {sheet_name} should have a recognised header, or every headerless row should have the same number of columns")
    column_count = row_lengths.pop()
    if column_count not in headerless_column_sets:
        expected_counts = ', '.join(str(count) for count in sorted(headerless_column_sets))
        raise ValueError(f"The {sheet_name} should have a recognised header or {expected_counts} headerless columns")
    return pd.DataFrame(rows, columns=headerless_column_sets[column_count])

def clean_string_columns(sheet):
    for column in sheet.columns:
        sheet[column] = sheet[column].fillna('').astype(str).str.strip()
    return sheet

def normalize_header_name(value):
    return str(value).strip().lower().strip('*`')

def is_provided(value):
    if pd.isna(value):
        return False
    return str(value).strip().lower() not in ['', 'na', 'none', 'null', '.']

def add_sample_identity(sample_sheet):
    input_file_ids = []
    sample_ids = []
    for _, row in sample_sheet.iterrows():
        input_file_id = input_file_id_from_path(row['reads'])
        sample_id = row.get('sample_id', 'NA')
        if not is_provided(sample_id):
            sample_id = infer_sample_id(input_file_id)
        input_file_ids.append(input_file_id)
        sample_ids.append(sample_id)
    sample_sheet['input_file_id'] = input_file_ids
    sample_sheet['sample_id'] = sample_ids
    return sample_sheet

def input_file_id_from_path(path):
    name = os.path.basename(str(path))
    for suffix in ['.fastq.gz', '.fq.gz', '.fasta.gz', '.fa.gz', '.fastq', '.fq', '.fasta', '.fa', '.bam']:
        if name.lower().endswith(suffix):
            return name[:-len(suffix)]
    return os.path.splitext(name)[0]

def infer_sample_id(input_file_id):
    return re.sub(r'(?i)(rep[0-9]+)[a-z]+(?=(_|$))', r'\1', input_file_id)

def validate_sample_sheet_genome_ids(sheet):
    if 'sample_sheet_genome_id' not in sheet.columns:
        return
    mismatched = sheet[
        sheet['sample_sheet_genome_id'].apply(is_provided) &
        (sheet['sample_sheet_genome_id'] != sheet['genome_id'])
    ]
    if not mismatched.empty:
        examples = [
            f"{row.reads}: sample sheet genome_id={row.sample_sheet_genome_id}, reference genome_id={row.genome_id}"
            for row in mismatched.head(5).itertuples(index=False)
        ]
        raise ValueError(
            "Sample-sheet genome_id values must match the genome_id resolved from the reference sheet. "
            f"Mismatches: {'; '.join(examples)}"
        )

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
    if has_missing_required(sample_sheet, ['reads', 'format', 'seq_source', 'map_key', 'species']):
        raise ValueError("The reads sample sheet has missing values in one of the required columns: reads, format, seq_source, map_key, species")
    
    # Check that the seq_source column is either 'pacbio' or 'nanopore'
    if not sample_sheet['seq_source'].isin(['pacbio', 'nanopore']).all():
        raise ValueError("The seq_source column in the reads sample sheet should have values either 'pacbio' or 'nanopore'")

    if not sample_sheet['format'].isin(['fasta', 'fastq', 'bam']).all():
        raise ValueError("The format column in the reads sample sheet should have values fasta, fastq, or bam")

    # read_id_map is optional provenance input for stable FASTA/FASTQ renaming.
    for _, row in sample_sheet.iterrows():
        read_id_map = row['read_id_map']
        if not is_provided(read_id_map):
            continue
        if row['format'] == 'bam':
            raise ValueError("read_id_map is only supported for FASTA/FASTQ rows; set it to NA for BAM rows")
        if not os.path.isfile(read_id_map):
            raise FileNotFoundError(f"The read_id_map file {read_id_map} does not exist")
        if not os.access(read_id_map, os.R_OK):
            raise PermissionError(f"The read_id_map file {read_id_map} is not readable")
    
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
    if has_missing_required(reference_sample_sheet, ['reference_file', 'genome_id', 'map_key', 'species']):
        raise ValueError("The reference sample sheet has missing values in one of the required columns: reference_file, genome_id, map_key, species")

    for reference_gtf in reference_sample_sheet['reference_gtf']:
        if not is_provided(reference_gtf):
            continue
        try:
            with open(reference_gtf) as f:
                pass
        except FileNotFoundError:
            raise FileNotFoundError(f"The reference GTF file {reference_gtf} does not exist")
        except Exception as e:
            raise Exception(f"The reference GTF file {reference_gtf} is corrupt")

    for _, row in reference_sample_sheet.iterrows():
        has_gtf = is_provided(row['reference_gtf'])
        has_polya = is_provided(row['polya_peak'])
        has_cage = is_provided(row['cage_peak'])
        if not has_gtf and (has_polya or has_cage):
            raise ValueError("The reference_gtf column is required when polya_peak or cage_peak is provided")
        if has_polya:
            validate_bed(row['polya_peak'], "polyA peak")
        if has_cage:
            validate_bed(row['cage_peak'], "CAGE peak")

    validate_sqanti3_reference_uniqueness(reference_sample_sheet)

def has_missing_required(sheet, columns):
    return any(not is_provided(value) for column in columns for value in sheet[column])

def validate_sqanti3_reference_uniqueness(reference_sample_sheet):
    sqanti_references = reference_sample_sheet[reference_sample_sheet['reference_gtf'].apply(is_provided)]
    for genome_id, rows in sqanti_references.groupby('genome_id', sort=True):
        evidence_tuples = rows[['reference_file', 'reference_gtf', 'polya_peak', 'cage_peak']].drop_duplicates()
        if len(evidence_tuples) <= 1:
            continue
        formatted = [
            "(" + ", ".join(str(value) for value in row) + ")"
            for row in evidence_tuples.itertuples(index=False, name=None)
        ]
        raise ValueError(
            f"Conflicting SQANTI3 reference rows found for genome_id '{genome_id}'. "
            "Rows with the same genome_id may use different map_key values, but "
            "reference_file, reference_gtf, polya_peak, and cage_peak must be identical. "
            f"Observed SQANTI3 input sets: {'; '.join(formatted)}"
        )

def validate_bed(path, label):
    try:
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt') as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 6:
                    raise ValueError(f"The {label} BED file {path} has fewer than 6 columns on line {line_number}")
                try:
                    int(fields[1])
                    int(fields[2])
                except ValueError:
                    raise ValueError(f"The {label} BED file {path} has non-integer start/end coordinates on line {line_number}")
                if fields[5] not in ['+', '-', '.']:
                    raise ValueError(f"The {label} BED file {path} has invalid strand '{fields[5]}' on line {line_number}")
                return
    except FileNotFoundError:
        raise FileNotFoundError(f"The {label} BED file {path} does not exist")
    except ValueError:
        raise
    except Exception as e:
        raise Exception(f"The {label} BED file {path} is corrupt")
    raise ValueError(f"The {label} BED file {path} has no data rows")

if __name__ == "__main__":
    main()
