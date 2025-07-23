#!/bin/bash
fasta="$1"
accession_id=${fasta/_fasta.gz/}
zcat "$fasta" | awk -v id="$accession_id" '/^>/ {print ">" id "_" ++i; next} {print}' > "$accession_id".reads.fasta
gzip -f "${accession_id}".reads.fasta
