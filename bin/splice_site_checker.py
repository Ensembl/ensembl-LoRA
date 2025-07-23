#!/usr/bin/env python 
import pysam
import logging
import argparse

def parse_bed12(bed12_file):
    with open(bed12_file, 'r') as f:
        for line in f:
            fields = line.strip().split()
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            trans_id = fields[3]
            score = fields[4]
            strand = fields[5]
            field_6 = fields[6]
            field_7 = fields[7]
            field_8 = fields[8]
            field_9 = fields[9]
            block_sizes = list(map(int, fields[10].split(',')))
            block_starts = list(map(int, fields[11].split(',')))
            yield chrom, start, end, trans_id, score, strand, field_6, field_7, field_8, field_9, block_sizes, block_starts

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))

def extract_splice_sites(chrom, start, strand, block_sizes, block_starts, fasta_file):
    splice_sites = []
    fasta = pysam.FastaFile(filename=fasta_file)
    for i in range(len(block_sizes) - 1):
        intron_start = start + block_starts[i] + block_sizes[i]
        intron_end = start + block_starts[i + 1]
        donor_site = fasta.fetch(chrom, intron_start, intron_start + 2).upper()
        acceptor_site = fasta.fetch(chrom, intron_end - 2, intron_end).upper()
        if strand == '-':
            donor_site = reverse_complement(donor_site)
            acceptor_site = reverse_complement(acceptor_site)
        splice_sites.append((donor_site, acceptor_site))
    return splice_sites

def check_canonical_splice_sites(splice_sites, strand):
    if not splice_sites:
        return False  # No splice sites means it's a single exon transcript
    if strand == '+':
        canonical_sites = [('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')]
    else:
        canonical_sites = [('AG', 'GT'), ('AG', 'GC'), ('AC', 'AT')]
    for donor, acceptor in splice_sites:
        if (donor, acceptor) in canonical_sites:
            return True
    return False

def process_bed12(bed12_file, fasta_file, prefix):
    for chrom, start, end, trans_id, score, strand, field_6, field_7, field_8, field_9, block_sizes, block_starts in parse_bed12(bed12_file):
        splice_sites = extract_splice_sites(chrom, start, strand, block_sizes, block_starts, fasta_file)
        if check_canonical_splice_sites(splice_sites, strand):
            logging.info(f"{trans_id} has canonical splice sites on {strand} strand: {splice_sites}")
            with open(f"{prefix}.fixed.bed", 'a', encoding='utf-8') as out_f:
                out_f.write(
                    f"{chrom}\t{start}\t{end}\t{trans_id}\t{score}\t{strand}\t"
                    f"{field_6}\t{field_7}\t{field_8}\t{field_9}\t"
                    f"{','.join(map(str, block_sizes))}\t"
                    f"{','.join(map(str, block_starts))}\n"
                )
        else:
            # Check the opposite strand
            opposite_strand = '+' if strand == '-' else '-'
            opposite_splice_sites = extract_splice_sites(chrom, start, opposite_strand, block_sizes, block_starts, fasta_file)
            if check_canonical_splice_sites(opposite_splice_sites, strand):
                logging.info(f"{trans_id} has Non-canonical splice sites on {strand} strand, but canonical on {opposite_strand} strand {splice_sites}. Swapping strand.")
                strand = opposite_strand
                splice_sites = opposite_splice_sites
                with open(f"{prefix}.fixed.bed", 'a', encoding='utf-8') as out_f:
                    out_f.write(f"{chrom}\t{start}\t{end}\t{trans_id}\t{score}\t{strand}\t{field_6}\t{field_7}\t{field_8}\t{field_9}\t{','.join(map(str, block_sizes))}\t{','.join(map(str, block_starts))}\n")
            elif not splice_sites:
                logging.info(f"{trans_id} is single exon. Removing transcript.")
                continue
            else:
                logging.info(f"{trans_id } has Non-canonical splice sites on both strands: {splice_sites}. Removing transcript.")
                continue  # Skip this transcript

def main():
    parser = argparse.ArgumentParser(description='Check splice sites in BED12 file.')
    parser.add_argument('bed12_file', help='Input BED12 file')
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('prefix', help='Prefix for bed12 output')
    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(filename=f"{args.prefix}.log", level=logging.INFO, format='%(message)s')

    # Process BED12 file
    process_bed12(args.bed12_file, args.fasta_file, args.prefix)

if __name__ == '__main__':
    main()