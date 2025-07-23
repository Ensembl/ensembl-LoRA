#!/usr/bin/env python 
import os
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="This script takes a SAM file and a bed file and adds strand information to the SAM file")
    parser.add_argument("--bam", help="BAM file", required=True)
    parser.add_argument("--prefix", help="Prefix for the output file", required=True)
    return parser.parse_args()

def main():
    args = get_args()

    # Extracting read and strand information from the merged bed file
    os.system('cat *.fixed.bed | sort -k1,1V -k2,2n | awk \'BEGIN { OFS="\t" } { print $4, $6 }\' > read_strand.txt')

    # Extract header from the SAM file
    os.system(f'samtools view -H {args.bam} > header.sam')
    
    # Filter out PG lines from the header
    os.system(f'grep -v "^@PG" header.sam > header.sam.tmp && mv header.sam.tmp header.sam')

    # Convert bam to sam file
    os.system(f'samtools view {args.bam} > temp.sam')

    # Read the strand information into a dictionary
    strand_info = {}
    with open("read_strand.txt", encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split("\t")
            strand_info[parts[0]] = parts[1]

    # Process the SAM file
    bam_file = 'temp.sam'
    with open(bam_file, encoding='utf-8') as f:
        with open('edited.sam', 'w', encoding='utf-8') as out_f:
            for line in f:
                if line.startswith("@"):
                    out_f.write(line)
                    continue
                parts = line.strip().split("\t")
                read_id = parts[0]
                if read_id in strand_info:
                    flag = int(parts[1])
                    if strand_info[read_id] == "+":
                        flag &= ~16  # Unset the 0x10 bit for positive strand
                    elif strand_info[read_id] == "-":
                        flag |= 16   # Set the 0x10 bit for negative strand
                    parts[1] = str(flag)
                    out_f.write("\t".join(parts) + "\n")

    # Add header to the processed SAM file
    os.system(f'cat header.sam edited.sam > {args.prefix}.canonical.sam')

    # Convert the processed SAM file to BAM
    os.system(f'samtools view -bS {args.prefix}.canonical.sam > {args.prefix}.canonical.bam')

    # Index the BAM file
    os.system(f'samtools index -c {args.prefix}.canonical.bam')

    # Remove temporary files
    os.system(f'rm temp.sam edited.sam read_strand.txt header.sam {args.prefix}.canonical.sam')

if __name__ == "__main__":
    main()
