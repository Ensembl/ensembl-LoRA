#!/usr/bin/env python
import argparse
import csv
import os
import re

import pandas as pd # type: ignore

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Convert tmerge output to GTF format")
    parser.add_argument("tmerge", help="tmerge output file")
    parser.add_argument("run_id", help="run id")
    parser.add_argument("output", help="gtf output file")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Assign the arguments to variables
    input = args.tmerge
    run_id = args.run_id
    output = args.output
    
    # Importing the tmerge output
    tmerge = pd.read_csv(f"{input}", sep="\t", header=None)
    
    # Converting tmerge output to tama gtf style
    # Extracting relevant columns and renaming them
    tmerge.columns = [
        "chromosome",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attributes",
    ]

    # Function to extract attribute values
    def extract_attributes(attr_string, attr_name):
        for attr in attr_string.split(";"):
            if attr.strip().startswith(attr_name):
                return attr.split('"')[1]
        return None

    # Extracting gene_id, transcript_id, exon_number, and exon_id
    tmerge["gene_id"] = tmerge["attributes"].apply(
        lambda x: extract_attributes(x, "gene_id")
    )
    tmerge["transcript_id"] = tmerge["attributes"].apply(
        lambda x: extract_attributes(x, "transcript_id")
    )
    tmerge["contains"] = tmerge["attributes"].apply(
        lambda x: extract_attributes(x, "contains")
    )

    # Replace all zeroes between _ and the first non-zero number
    tmerge["gene_id"] = tmerge["gene_id"].apply(
        lambda x: re.sub(r"_(0+)([1-9])", r"_\2", x)
    )
    tmerge["transcript_id"] = tmerge["transcript_id"].apply(
        lambda x: re.sub(r"_(0+)([1-9])", r"_\2", x)
    )

    # Adding exon_number and exon_id manually
    tmerge["exon_number"] = tmerge.groupby("transcript_id").cumcount() + 1
    tmerge["exon_id"] = tmerge.apply(
        lambda row: f"{row['transcript_id']}.{row['exon_number']}", axis=1
    )
    # Dropping the original attributes column
    tmerge = tmerge.drop(columns=["attributes"])

    # Editing gene, transcript and exon ids
    tmerge["gene_id"] = tmerge["gene_id"].apply(
        lambda x: x.replace("chr", f"{run_id}_", 1)
    )
    tmerge["transcript_id"] = tmerge["transcript_id"].apply(
        lambda x: x.replace("chr", f"{run_id}_", 1)
    )
    tmerge["exon_id"] = tmerge["exon_id"].apply(
        lambda x: x.replace("chr", f"{run_id}_", 1)
    )

    # Extracting the transcript number from transcript_id
    tmerge["transcript_number"] = tmerge["transcript_id"].apply(
        lambda x: int(x.split("_")[-1])
    )

    # Sorting based on the transcript number
    tmerge = tmerge.sort_values(by="transcript_number")

    # Removing the transcript_number column
    tmerge = tmerge.drop(columns=["transcript_number"])

    # Creating contains file
    tmerge[["chromosome", "transcript_id", "contains"]].drop_duplicates().to_csv(
        f"{output}.contains.tsv", sep="\t", header=False, index=False
    )

    # Merging all attributes
    tmerge["attributes"] = tmerge.apply(
        lambda row: f'gene_id "{row["gene_id"]}"; transcript_id "{row["transcript_id"]}"; exon_number "{row["exon_number"]}"; exon_id "{row["exon_id"]}";',
        axis=1,
    )

    # Removing columns
    tmerge = tmerge[
        [
            "chromosome",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attributes",
        ]
    ]

    # Saving the modified tmerge output
    tmerge.to_csv(
        f"{output}.tama.gtf",
        sep="\t",
        header=False,
        index=False,
        quoting=csv.QUOTE_NONE,
    )

if __name__ == '__main__':
    main()