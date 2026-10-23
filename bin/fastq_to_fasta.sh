#!/usr/bin/env bash 
fastq="$1"
output="$2"
seqtk seq -A "$fastq" > "$output"
gzip "$output"
if [ $? -eq 0 ]; then
    rm "$fastq" $(readlink -f "$fastq")
fi