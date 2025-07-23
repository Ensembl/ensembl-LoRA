#!/usr/bin/env python
import pandas as pd
import argparse

# Calculate intron chains
def calculate_intron_chains(group):
    strand = group['strand'].values[0]
    if strand == "+" or strand == 1:
        intron_starts = group['end'].values[:-1] + 1
        intron_ends = group['start'].values[1:] - 1
        return ",".join(f"{start}-{end}" for start, end in zip(intron_starts, intron_ends))
    elif strand == "-" or strand == -1:
        intron_starts = group['end'].values[1:] + 1
        intron_ends = group['start'].values[:-1] - 1
        return ",".join(f"{start}-{end}" for start, end in zip(intron_starts, intron_ends))

def main():
    parser = argparse.ArgumentParser(description='Extract intron and exon information from proto-transcript files')
    parser.add_argument('--transcripts', '-t', type=str, help='GTF/GFF3 containing the proto-transcripts')
    parser.add_argument('--prefix', '-p', type=str, help="Output file prefix")
    args = parser.parse_args()
    
    # Load the GTF/GFF3 file
    proto_transcripts = []
    with open(args.transcripts, "r") as transcript_file:
        for exon in transcript_file:
            columns = exon.strip().split("\t")
            transcript_id = columns[-1].split(";")[1].split(" ")[-1].replace('"', '')
            processed_line = [columns[0], columns[3], columns[4], columns[6], transcript_id]
            proto_transcripts.append(processed_line)
    proto_transcripts = pd.DataFrame(proto_transcripts, columns=["seqname", "start", "end", "strand", "transcript_id"])

    # Preparing start and end columns
    proto_transcripts['start'] = proto_transcripts['start'].astype(int)
    proto_transcripts['end'] = proto_transcripts['end'].astype(int)

    # Split by strand the GTF file
    input_positive = proto_transcripts[((proto_transcripts['strand'] == 1) | (proto_transcripts['strand'] == "+"))]
    input_negative = proto_transcripts[((proto_transcripts['strand'] == -1) | (proto_transcripts['strand'] == "-"))]

    # Sort by transcript and start position
    input_positive = input_positive.sort_values(by=['transcript_id', 'start'])
    input_negative = input_negative.sort_values(by=['transcript_id', 'start'], ascending = False)

    # Combine the two dataframes
    input_combined = pd.concat([input_positive, input_negative])

    # Calculate intron chains
    intron_chains = input_combined.groupby(['seqname', 'transcript_id']).apply(calculate_intron_chains, include_groups=False).reset_index(name='intron_chain')
    del input_positive, input_negative, input_combined
    intron_chains.loc[:, "intron_chain"] = "(" + intron_chains["seqname"].astype(str) + ")" + intron_chains["intron_chain"].astype(str)
    single_exon_transcripts = intron_chains[intron_chains['intron_chain'] == intron_chains['seqname'].apply(lambda x: f'({x})')].reset_index(drop=True)
    multi_exon_transcripts = intron_chains[intron_chains['intron_chain'] != intron_chains['seqname'].apply(lambda x: f'({x})')].reset_index(drop=True)
    # Adding strand information to multiple exon transcripts
    multi_exon_transcripts = multi_exon_transcripts.merge(proto_transcripts[["seqname", "transcript_id", "strand"]].drop_duplicates(), on=["seqname", "transcript_id"], how="left")
    multi_exon_transcripts["intron_chain"] = multi_exon_transcripts["strand"].apply(lambda x: f"({x})") + multi_exon_transcripts["intron_chain"]
    multi_exon_transcripts = multi_exon_transcripts[['seqname', 'transcript_id', 'intron_chain']]

    # Extracting exon information for single exon transcripts
    single_exons = proto_transcripts[proto_transcripts['transcript_id'].isin(single_exon_transcripts['transcript_id'])]

    # Splitting single exon transcripts by strand
    single_exons_positive = single_exons[((single_exons['strand'] == 1) | (single_exons['strand'] == "+"))]
    single_exons_negative = single_exons[((single_exons['strand'] == -1) | (single_exons['strand'] == "-"))]

    # Preparing exon information
    single_exons_positive = single_exons_positive.copy()
    single_exons_negative = single_exons_negative.copy()
    single_exons_positive.loc[:, "exon_info"] = "(" + single_exons_positive["strand"].astype(str) + ")" + "(" + single_exons_positive["seqname"].astype(str) + ")" + single_exons_positive['start'].astype(str) + "-" + single_exons_positive['end'].astype(str)
    single_exons_negative.loc[:, "exon_info"] = "(" + single_exons_negative["strand"].astype(str) + ")" + "(" + single_exons_negative["seqname"].astype(str) + ")" + single_exons_negative['end'].astype(str) + "-" + single_exons_negative['start'].astype(str)

    # Combine the two dataframes
    single_exons = pd.concat([single_exons_positive, single_exons_negative])
    single_exons = single_exons[["seqname", "transcript_id", "exon_info"]]

    # Free up memory by removing large dataframes
    del proto_transcripts, single_exons_positive, single_exons_negative, intron_chains

    # Save intron/exon chains to TSV
    multi_exon_transcripts.to_csv(f'{args.prefix}_intron_chains.tsv', sep="\t", index=False, header=False)
    del multi_exon_transcripts
    single_exons.to_csv(f'{args.prefix}_single_exons.tsv', sep="\t", index=False, header=False)
    del single_exons

if __name__ == "__main__":
    main()