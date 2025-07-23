#!/bin/bash
fastq="$1"
accession_id=${fastq%_fastq.gz}
zcat "$fastq" | awk '{if(NR%4==1) {print "@'"$accession_id"'_" NR/4 + 0.75} else {print $0}}' > "${accession_id}.reads.fastq"
gzip -f "${accession_id}.reads.fastq"
