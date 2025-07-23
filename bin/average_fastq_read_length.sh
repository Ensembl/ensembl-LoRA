#!/usr/bin/env bash 
fastq="$1"
awk_output=$(zcat "${fastq}" | awk 'NR%40000==2 {count++; sum+=length} END {printf "%.0f", sum/count}')

if (( awk_output > 200 )); then
    cat ${fastq} > ${fastq%.fastq.gz}_fastq.gz
else
    cat ${fastq} > ${fastq%.fastq.gz}.failed
fi
