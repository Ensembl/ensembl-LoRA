#!/usr/bin/env perl

#Determines the correct strand using canonical splice sites and polyA sequence
#Input: BED12 file transformed from a BAM file
#Output: BED12 file with strand values fixed 

use strict;
use Getopt::Long;
use Bio::DB::HTS::Faidx;
$|=1;

my $infile;
my $outfile;
my $genome_fasta;
my $verbose = 0;

&GetOptions(
            'infile=s'  => \$infile,
            'outfile=s' => \$outfile,
            'fasta=s'   => \$genome_fasta,
            'verbose!'  => \$verbose
           );

my $index = Bio::DB::HTS::Faidx->new($genome_fasta);

open (OUT, ">$outfile") or die "Can't open file $outfile\n";
open (IN, $infile) or die "Can't open file $infile\n";
while (<IN>){
  #chr13	16327520	16440742	6214ab95-faea-464e-ab04-8612ba37d8e2	0	-	16327520	16440742	55,0,0	4	251,543,334,196	0,19363,93590,113026
  next if /^#/;
  chomp;
  my @fs = split(/\t/);
  my $chr = $fs[0];
  my $start = $fs[1];
  my $end = $fs[2];
  my $name = $fs[3];
  my $strand = $fs[5];
  my @exon_lengths = split(/,/, $fs[10]);
  my @exon_starts = split(/,/, $fs[11]);
  if (scalar @exon_starts > 1){
    #Look at the canonical splice sites
    my ($right_strand, $wrong_strand, $unknown) = (0, 0, 0);
    print OUT "#$name ($strand)  " if $verbose;
    for (my $i =0; $i < scalar(@exon_starts) - 1; $i++){
      my $spl1 = $start + $exon_starts[$i] + $exon_lengths[$i];
      my $spl2 = $start + $exon_starts[$i+1] + 1;
      my $seq1 = $index->get_sequence_no_length($chr.":".($spl1 + 1)."-".($spl1 + 2));
      my $seq2 = $index->get_sequence_no_length($chr.":".($spl2 - 2)."-".($spl2 - 1));
      my $splice_site_seq = $seq1."..".$seq2; 
      print OUT $splice_site_seq."  " if $verbose;
      my $seq_type = "other";
      if ($splice_site_seq eq "GT..AG" or $splice_site_seq eq "GC..AG" or $splice_site_seq eq "AT..AC"){
        $seq_type = "canonical";
      }
      elsif ($splice_site_seq eq "CT..AC" or $splice_site_seq eq "CT..GC" or $splice_site_seq eq "GT..AT"){
        $seq_type = "rev_comp_canonical";
      }
      else{
        $seq_type = "non_canonical";  
      }
	  
      if (($strand eq "+" and $seq_type eq "canonical") or ($strand eq "-" and $seq_type eq "rev_comp_canonical")){
        $right_strand++; 
      }
      elsif (($strand eq "-" and $seq_type eq "canonical") or ($strand eq "+" and $seq_type eq "rev_comp_canonical")){
        $wrong_strand++;
      }
      else{
        $unknown++;  
      }
    } 
 
    if (!($right_strand) and $wrong_strand >= $unknown){
      #print line with corrected strand
      print OUT "CORRECTED\n" if $verbose;
      $fs[5] = ($strand eq "+" ? "-" : "+");
      print OUT join("\t", @fs)."\n";
    }
    else{
      #print the original line
      print OUT "OK\n" if $verbose;
      print OUT $_."\n";
    }		
  }
  else{
    #TO DO: look at the polyA / polyT tail
    print OUT "#SINGLE_EXON\n" if $verbose;
    print OUT $_."\n"; 
  }	
}
close (IN);
close (OUT);


