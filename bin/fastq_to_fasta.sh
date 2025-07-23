#!/usr/bin/env bash 
fastq="$1"
output="$2"
seqtk seq -A "$fastq" > "$output"
gzip "$output"
