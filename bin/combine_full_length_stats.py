#!/usr/bin/env python 

import pandas as pd
import os
import argparse


def main():
    args = get_args()
    files_to_import = args.files
    imported_files = []
    for file in files_to_import:
        imported_files.append(pd.read_csv(file, sep='\t'))
    df = pd.concat(imported_files)
    df.to_csv(args.prefix + '_full_length_stats.tsv', sep='\t', index=False)

def get_args():
    parser = argparse.ArgumentParser(description='Combine full length stats files')
    parser.add_argument('-f', '--files', type=str, nargs='*', help='Files to combine')
    parser.add_argument('-p', '--prefix', type=str, help='Prefix for output file')
    return parser.parse_args()

if __name__ == "__main__":
    main()