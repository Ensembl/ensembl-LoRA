#!/usr/bin/env python
import re
import argparse
import time
import os

# Create the parser
parser = argparse.ArgumentParser(description='Produce usable Nextflow progress report from slurm output file')

# Add the arguments
parser.add_argument('-s', '--shard', default='', nargs='+', required=True, help='Input shard files')
parser.add_argument('-b', '--bam', default='', required=True, help='Input Bam file for name replacement')

# Execute the parse_args() method
argument = parser.parse_args()

# Get the shard files
shard_files = argument.shard

# Get the bam file
bam_file = argument.bam

# Get the name of the bam file
bam_name = os.path.basename(bam_file).split('.')[0]
bam_name = f'{bam_name}_'

# Get the name of the shard files
output_shards = [f'{bam_name}{os.path.basename(shard_file).split(".")[0].split("_")[1]}.picard.bam' for shard_file in shard_files]

# Rename the shard files
for shard_file, output_shard in zip(shard_files, output_shards):
    os.rename(shard_file, output_shard)