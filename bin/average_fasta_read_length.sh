#!/usr/bin/env bash 
fasta="$1"
awk_output=$(zcat "${fasta}" | awk '/^>/ {if (seqlen){count++; sum+=seqlen}; seqlen=0; next} {seqlen += length} END {if (seqlen){count++; sum+=seqlen}; printf "%.0f", sum/count}')

if (( awk_output > 200 )); then
    cat ${fasta} > ${fasta%.fasta.gz}_fasta.gz
else
    cat ${fasta} > ${fasta%.fasta.gz}.failed
fi