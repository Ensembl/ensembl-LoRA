#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
no warnings 'recursion';

my $message_text  = "Error\n";
my $exit_status   = 2;          ## The exit status to use
my $verbose_level = 99;          ## The verbose level to use
my $filehandle    = \*STDERR;   ## The filehandle to write to
my $sections = "NAME|SYNOPSIS|DESCRIPTION";


=head1 NAME

tmerge

=head1 SYNOPSIS

Merge transcriptome read-to-genome alignments into non-redundant transcript models.

C<tmerge> compares transcript structures (or read-to-genome alignments) present in the input and attempts to reduce transcript redundancy, I<i.e.>, merge compatible input transcripts into non-redundant transcript models. The program treats spliced and monoexonic reads separately (I<i.e.>, those are never merged together).

C<tmerge> is fast and can typically process several millions of aligned long reads in a few minutes.

=begin HTML

<p><img src="http://public-docs.crg.es/rguigo/CLS/img/tmerge1.50dpi.png" alt="tmerge sketch" /></p>

=end HTML

See DESCRIPTION below for more details.

B<Usage example>:

C<< tmerge --tmPrefix <custom transcript_id prefix string for output GTF> <input GTF file> > <output file> >>


=head2 INPUT

GTF file of read-to-genome alignments, sorted by chromosome and start position.

Only C<exon> records are considered.
Read alignments need to be uniquely identified with the C<transcript_id> GTF attribute. C<transcript_id> is the only mandatory GTF attribute in input records.

=head2 OPTIONS

=over

=item * C<tmPrefix> (string) = Prefix string for C<transcript_id> identifiers in the output

B<Default>: '' (empty string)

By default, output C<transcript_id>s consist in arbitrary "C<TM_XXXXXXXXXXXX>" strings. If C<tmPrefix> is set, its value will prefix all C<transcript_id> strings in the GTF output.

=item * C<minReadSupport> (integer) = minimum number of times a read alignment (as defined by its exon/intron structure) needs to be present in the input. In other words, when building a transcript model, only the reads fulfilling the following conditions are considered:

For B<spliced transcripts>, at least C<minReadSupport> input reads must share a given intron chain and 5' + 3' ends (+/- C<endFuzz> bases, see below).

For B<mono-exonic transcripts>, at least C<minReadSupport> input reads must share their 5' + 3' ends (+/- C<endFuzz> bases, see below). In other words, when C<endFuzz> C< = 0> (the default), only monoexonic reads with identical genome coordinates are merged.

B<Default>: 1

=item * C<endFuzz> (positive integer) = Tolerated fuzziness of 5' and 3' ends for two reads to be considered equivalent when calculating read support (see C<minReadSupport> option above)

B<Default>: 0 (i.e., no fuzziness allowed)


=item * C<exonOverhangTolerance> (positive integer) = maximum number of nucleotides of terminal exon overhang allowed within an intron of another transcript during the merging of input reads. See explanation in "DESCRIPTION" below.

B<Default>: 0 (i.e., no exon overhang allowed)

=back

=head2 OUTPUT

C<tmerge> outputs non-redundant transcript models (B<TMs>) in GTF format. Each TM entry is uniquely identified by its (arbitrary) C<transcript_id> attribute.

The C<gene_id> attribute has the same value as C<transcript_id> by convention; it is therefore meaningless.

The following extra GTF attributes are present in the 9th field, in order:

=over

=item * C<contains> (string): comma-separated list of input reads (C<transcript_id>s) contained in the TM, sorted by descending genomic size.

=item * C<contains_count> (integer): number of input reads contained in the TM.

=item * C<3p_dists_to_3p> (string): comma-separated list of the distances (always positive, in bases on mature RNA, i.e. ignoring introns) of the TM's 3' end to each of the input reads 3' ends it C<contains>. The list's order follows that of C<contains>.

=item * C<5p_dists_to_5p> (string): comma-separated list of the distances (always positive, in bases on mature RNA, i.e. ignoring introns) of the TM's 5' end to each of the input reads 5' ends it C<contains>. The list's order follows that of C<contains>.

=item * C<flrpm> (float): TM's expression quantification in "Full-Length Reads per Million". This corresponds to C<longest_FL_supporters_count> divided by the number of reads (i.e., C<transcript_id>'s) present in the input.

=item * C<longest> (string): comma-separated list of the longest read(s) (C<transcript_id>s) contained in the TM. This list contains more that one item only in case of length ties. Note that the reads reported do not necessarily cover the entire length of the resulting TM.

=item * C<longest_FL_supporters> (string): comma-separated list of input reads that support C<longest> over C<longest>'s full-length (+/- C<endFuzz>).

=item * C<longest_FL_supporters_count> (integer): number of input reads that support C<longest> over C<longest>'s full-length (+/- C<endFuzz>).

=item * C<mature_RNA_length> (integer): the mature RNA length of the TM (i.e., the sum of the lengths of all its exons)

=item * C<meta_3p_dists_to_5p> (string): comma-separated list of the distances (comprised between 0 and 1, on mature RNA, i.e. ignoring introns) of the TM's B<5' end> to each of the input reads 3' ends it C<contains>, normalized over the TM's mature RNA length. The list's order follows that of C<contains>.

=item * C<meta_5p_dists_to_5p> (string): comma-separated list of the distances (comprised between 0 and 1, on mature RNA, i.e. ignoring introns) of the TM's 5' end to each of the input reads 5' ends it C<contains>, normalized over the TM's mature RNA length. The list's order follows that of C<contains>.

=item * C<rpm> (float): TM's expression quantification in "Reads per Million". This corresponds to C<contains_count> divided by the number of reads (i.e, C<transcript_id>'s) present in the input.

=item * C<spliced> (boolean): specifies if the TM is spliced (1) or monoexonic (0).

=back

=head1 DESCRIPTION

C<tmerge> reduces redundancy in a set of transcriptome read-to-genome alignments. It does so by looking for reads with I<B<compatible>> aligned structures in the input, and merging those into I<B<Transcript Models>> (B<TMs>).

=begin HTML

<p><img src="http://public-docs.crg.es/rguigo/CLS/img/tmerge1.50dpi.png" alt="tmerge sketch" /></p>

=end HTML

Pairwise B<compatibility> between aligned structures is evaluated using the following rules:

=over

=item * If both structures are B<spliced>, they are deemed compatible if:

=over

=item * 1. at least one of their exons overlap on the same genomic strand,

=item * 2. either their intron chains are equal, or one is an exact subset of the other,

and

=item * 3. there is no overlap between an exon of one structure and an intron of the other.

=back

Condition (2) means that C<tmerge> will never artificially extend intron chains:

=begin HTML

<p><img src="http://public-docs.crg.es/rguigo/CLS/img/tmerge2.50dpi.png" alt="tmerge non-merge case" /></p>

=end HTML

=item * If both structures are B<monoexonic>, they are considered compatible if they overlap by at least 1 nucleotide on the same genomic strand.

=item * If one structure is B<spliced> and the other B<monoexonic>, they are not merged.

=back


All pairs of compatible structures are then merged recursively into the longest possible TM.

=head2 C<exonOverhangTolerance> option and splice sites

Setting this option to a positive integer can correct mismapped splice junctions that sometimes occur when aligning very short, error-rich terminal read exons:


=begin HTML

<p><img src="http://public-docs.crg.es/rguigo/CLS/img/tmerge_FalseExonOverhang.50dpi.png" alt="tmerge FalseExonOverhang sketch" /></p>

=end HTML

The setting works as explained below:

=begin HTML

<p><img src="http://public-docs.crg.es/rguigo/CLS/img/tmerge_exonOverhangTolerance.50dpi.png" alt="tmerge exonOverhangTolerance sketch" /></p>

=end HTML


=head1 AUTHOR

Julien Lagarde, CRG, Barcelona, contact julienlag@gmail.com

=cut



my $tmPrefix='';
my $minReadSupport=1;
my $exonOverhangTolerance=0;
my $transcriptEndFuzziness=0;
my $debug='';
GetOptions ('tmPrefix=s' => \$tmPrefix,
            'minReadSupport=i' => \$minReadSupport,
            'exonOverhangTolerance=i' => \$exonOverhangTolerance,
						'endFuzz=i' => \$transcriptEndFuzziness,
            'debug' => \$debug,
            )
or pod2usage( { -message => "Error in command line arguments",
        		  -exitval => $exit_status  ,
            		-verbose => $verbose_level,
               -output  => $filehandle } );

unless (defined $ARGV[0]){
	pod2usage( { -message => "Error in command line arguments: no input provided.",
        		  -exitval => $exit_status  ,
            		-verbose => $verbose_level,
               -output  => $filehandle } );
}

die "ERROR: minReadSupport value must be > 0, can't continue.\n" if $minReadSupport<1;
die "ERROR: exonOverhangTolerance value must be >= 0, can't continue.\n" if $exonOverhangTolerance<0;
die "ERROR: endFuzz value must be >= 0, can't continue.\n" if $transcriptEndFuzziness<0;

print STDERR "######################################################################################################################
######################################################################################################################
##
";

print STDERR "##    Parameters information:\n";
print STDERR "##\n";
print STDERR "##    - TM transcript_id output prefix (--tmPrefix): $tmPrefix\n";
print STDERR "##    - Minimum read support required per merged transcript model (--minReadSupport): $minReadSupport\n";
print STDERR "##    - Exon overhang tolerance when merging (--exonOverhangTolerance): $exonOverhangTolerance bases\n";
print STDERR "##    - End fuzziness tolerance when calculating support (--endFuzz): $transcriptEndFuzziness bases\n";
print STDERR "##\n";
print STDERR "######################################################################################################################
######################################################################################################################
";


my $sortedGff;
if(defined($ARGV[0])){
	$sortedGff=$ARGV[0];
}
else{
	die "ERROR: Need input GTF file name as argument.\n";
}

open GFF, "$sortedGff" or die "ERROR: ".$!;

my %transcript_to_transcript=();
my %transcript_exons=();
my %transcript_chr=();
my %transcript_strand=();
my $previous_start=-1;
my $previous_chr='Caravaggio';
my $previous_transcript='Picasso';
my $superExonStart=-1;
my $superExonStop=-1;
my %transcript_id_index=();

print STDERR "Parsing GTF input...\n";
my $nr_exons=0;
my $read_count=0;
my %transcript_seen=();
my %transcript_index=();
my %transcript_rev_index=();
while (<GFF>){
	next if ($_=~/^#/);
	next unless ($_=~/\texon\t/);
	if ($_=~/^(\S+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t\S+\t(\S+)\t\S+\t.*transcript_id "(\S+)?";/){
		if ($3 eq "exon"){
			$nr_exons++;
			my $GTFtranscript_id=$7;
			my $chr=$1;
			my $start=$4;
			my $stop=$5;
			my $str=$6;
			my $strand;
			if($str eq '+'){
				$strand=1;
			}
			elsif($str eq '-'){
				$strand=-1;
			}
			elsif($str eq '.'){
				$strand=0;
			}
			else{
				die "ERROR: Unrecognized strand value '$str' at line $.\n";
			}
      die "ERROR: Corrupt GTF. Start coordinate cannot be greater than stop coordinate at line $. . Can't continue.\n" if($start > $stop);
			unless(exists $transcript_seen{$GTFtranscript_id}){
				$read_count++;
				$transcript_seen{$GTFtranscript_id}=undef;
				$transcript_index{$read_count}=$GTFtranscript_id;
				$transcript_rev_index{$GTFtranscript_id}=$read_count;
			}

			my $transcript_id=$transcript_rev_index{$GTFtranscript_id};
			#Check for sorted input:
			die "ERROR: Unsorted GTF input (line $.). Must be sorted by chr, then start, then stop. Can't continue.\n" if ($chr eq $previous_chr && $start < $previous_start);
      if(exists $transcript_strand{$transcript_id} && $transcript_strand{$transcript_id} != $strand){
        die "ERROR: Inconsistent strand for transcript $GTFtranscript_id in input file. Can't continue.\n";
      }
			$transcript_strand{$transcript_id}=$strand;
      if (exists $transcript_chr{$transcript_id} && $transcript_chr{$transcript_id} ne $chr){
        die "ERROR: Inconsistent chr for transcript $GTFtranscript_id in input file. Can't continue.\n";

      }
			$transcript_chr{$transcript_id}=$chr;
			my @exon=($chr, $start, $stop, $strand);
			push(@{$transcript_exons{$transcript_id}}, \@exon);
			$transcript_to_transcript{$transcript_id}{$transcript_id}=undef;

			### BEGIN this should be a subroutine but repeated subroutine calls are too expensive in perl
			### my $overlapWithPrevious=overlap($superExonStart, $superExonStop, $start, $stop);
			my $start1=$superExonStart;
			my $stop1=$superExonStop;
			my $start2=$start;
			my $stop2=$stop;
			my $overlap;
			my $start2minusstop1=$start2-$stop1;
			my $start2minusstart1=$start2-$start1;
			my $stop2minusstart1=$stop2-$start1;
			my $stop2minusstop1=$stop2-$stop1;
			if( ( $stop2minusstart1>=0  && $stop2minusstop1 <=0 ) || ($start2minusstart1 >=0 && $start2minusstop1 <=0) || ($start2minusstart1 <= 0 && $stop2minusstop1 >= 0)){
					$overlap=1;
			}
			else{
				if($stop2minusstart1<0){
					$overlap=-1
				}
				elsif($start2minusstop1>0){
					$overlap=-2
				}
			}
			### # 1  : overlap
	            # -1 : 2 upstream of 1
	            # -2 : 2 downstream of 1
			### END this should be a subroutine but repeated subroutine calls are too expensive in perl




			if ($chr eq $previous_chr && $overlap == 1){
				$transcript_to_transcript{$transcript_id}{$previous_transcript}=undef;
				$transcript_to_transcript{$previous_transcript}{$transcript_id}=undef;
				if($stop > $superExonStop){
					$superExonStop = $stop;
				}
			}

			else{
				#re-initialize superExon
				$previous_chr=$chr;
				$superExonStart=$start;
				$superExonStop=$stop;
			}
			$previous_transcript=$transcript_id;
			$previous_start=$start;
		}
	}
	else{
		die "ERROR: line $.: malformed GTF record. Can't continue.\n";

	}

}
close GFF;
print STDERR "Done. Found $nr_exons exons and $read_count transcripts.\n";
my $million_read_count=1000000/($read_count+1);
%transcript_seen=();
%transcript_rev_index=();

print STDERR "Building contigs (sets of overlapping transcripts)...\n";
#  we build contigs to reduce the search space when looking for compatible transcript structures
my $locusNumber = 0;
my %transcript_id_to_locus_id=();
foreach my $tr1 (keys %transcript_to_transcript){
	buildContig($tr1, \%transcript_to_transcript, \%transcript_id_to_locus_id, $locusNumber);
	#generate UUID for contig (this is to write temp files):
	$locusNumber++;
}

%transcript_to_transcript=(); #free up some memory

my %contig_to_transcripts=();
foreach my $tr (keys %transcript_id_to_locus_id){
	my $contig=$transcript_id_to_locus_id{$tr};
	push(@{$contig_to_transcripts{$contig}}, $tr);
}
%transcript_id_to_locus_id=();
my $nrContigs=scalar(keys %contig_to_transcripts);
print STDERR "Done. Built $nrContigs contigs.\n";

print STDERR "Comparing transcript structures...\n";
my %transcript_introns=();
my %container_original_ends=(); #contains non-adjusted end coordinates for all containers

CONTIGS: foreach my $contig (keys %contig_to_transcripts){
	my %container_to_transcripts=(); #'value' is a subset of , or equal to, 'key''s intron chain.
  my %container_to_supporting_transcripts=(); # same as %container_to_transcripts but more stringent ('key' contains 'value's which consist only in equal transcripts +/- exonOverhangTolerance)
	my %transcript_to_container=();
	my @list1=@{$contig_to_transcripts{$contig}};
	#build introns within contig:
	my @list1Mono=();
	my @list1Spliced=();
	my %supportCount=();
	# build sets of introns for each transcript, populate @list1Mono and/or @list1Spliced accordingly:
	foreach my $tr (@list1){
		$supportCount{$tr}=1;
		for (my $i=0; $i< $#{$transcript_exons{$tr}};$i++){
			my $intronChr=${$transcript_exons{$tr}}[$i][0];
			my $intronStrand=${$transcript_exons{$tr}}[$i][3];
			my $intronStart=${$transcript_exons{$tr}}[$i][2]+1;
			my $intronStop=${$transcript_exons{$tr}}[$i+1][1]-1;
			my @intron=($intronChr, $intronStart, $intronStop, $intronStrand);
			push(@{$transcript_introns{$tr}}, \@intron);

		}
		if($#{$transcript_exons{$tr}} == 0){
			push (@list1Mono, $tr)
		}
		else{
			push(@list1Spliced, $tr);
		}
	}
	@list1=();

	@list1Spliced= sort ({ ${$transcript_introns{$a}}[0][1] <=> ${$transcript_introns{$b}}[0][1] or ${$transcript_exons{$a}}[0][1] <=> ${$transcript_exons{$b}}[0][1] or ${$transcript_exons{$b}}[-1][2] <=> ${$transcript_exons{$a}}[-1][2] or $a cmp $b } @list1Spliced); #sort spliced transcripts by position of first intron. second and third comparisons are necessary so the script is deterministic (otherwise ties are handled randomly)

	# compute read support for each input spliced transcript/read:
	print STDERR "Calculating spliced read support...\n" if $debug;
	NEXTTR1: for (my $k=0; $k<=$#list1Spliced; $k++){
		my $tr1=$list1Spliced[$k];
		print STDERR "DEBUG: tr1 $transcript_index{$tr1}\n" if $debug;
		if(exists ($transcript_to_container{$tr1})){
			#$tr1=$transcript_to_container{$tr1};
			print STDERR "DEBUG: tr1 IS CONTAINED in $transcript_index{$transcript_to_container{$tr1}}\n" if $debug;
			next;
		}
    @{$container_original_ends{$tr1}}=(${$transcript_exons{$tr1}}[0][1], ${$transcript_exons{$tr1}}[-1][2]);

		NEXTTR2: for (my $l=$k+1; $l<=$#list1Spliced; $l++){
			my $tr2=$list1Spliced[$l];
			print STDERR "DEBUG:  tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2}\n" if $debug;
			if(exists ($transcript_to_container{$tr2})){
			#$tr2=$transcript_to_container{$tr2};
				next;
			}
			if(${$transcript_exons{$tr1}}[-1][2] - ${$transcript_exons{$tr2}}[0][1]>0){
				if( ( $transcript_strand{$tr1} == $transcript_strand{$tr2})
				   &&
				   ($#{$transcript_introns{$tr2}} == $#{$transcript_introns{$tr1}})
				   &&
				   (${$transcript_exons{$tr2}}[0][1] >= ${$transcript_exons{$tr1}}[0][1] - $transcriptEndFuzziness
				   &&
				   ${$transcript_exons{$tr2}}[0][1] <= ${$transcript_exons{$tr1}}[0][1] + $transcriptEndFuzziness)
				   &&
				   (${$transcript_exons{$tr2}}[-1][2] >= ${$transcript_exons{$tr1}}[-1][2] - $transcriptEndFuzziness
				    &&
				    ${$transcript_exons{$tr2}}[-1][2] <= ${$transcript_exons{$tr1}}[-1][2] + $transcriptEndFuzziness) ){
					print STDERR "DEBUG:  tr1 $transcript_index{$tr1} overlaps tr2 $transcript_index{$tr2}\n" if $debug;

					for (my $j=0; $j <= $#{$transcript_introns{$tr2}}; $j++){
						for (my $i=$j; $i<= $#{$transcript_introns{$tr1}}; $i++){
							if(${${$transcript_introns{$tr2}}[$j]}[1] == ${${$transcript_introns{$tr1}}[$i]}[1] && ${${$transcript_introns{$tr2}}[$j]}[2] == ${${$transcript_introns{$tr1}}[$i]}[2]){
								print STDERR "DEBUG:    tr1 intron $i = tr2 intron $j\n" if $debug;
								last;
							}
							else{
								print STDERR "DEBUG:    tr1 intron $i != tr2 intron $j\n" if $debug;
								next NEXTTR2;
							}
						}
					}
					$container_to_transcripts{$tr1}{$tr2}=undef ;
					$transcript_to_container{$tr2}=$tr1;
					adjustContainerEnds($tr1,$tr2);
					if (exists $container_to_transcripts{$tr2}){
						foreach my $tr3 (keys %{$container_to_transcripts{$tr2}}){
							$transcript_to_container{$tr3}=$tr1;
							$container_to_transcripts{$tr1}{$tr3}=undef ;
						}
						delete($container_to_transcripts{$tr2});
					}


					print STDERR "DEBUG:     tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2} FULLMATCH 2\n" if $debug;
					$supportCount{$tr1}++;
					$container_to_supporting_transcripts{$tr1}{$tr2}=1;
				}
			}
			else{ #tr2 is downstream of tr1. Since trs are sorted by start position, we can safely skip to next tr1
		 		print STDERR "DEBUG:   tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2} INCOMP 10\n" if $debug;
		 		last;
			}
		}
	}

	print STDERR "Merging spliced reads...\n" if $debug;

	NEXTTR1: for (my $k=0; $k<=$#list1Spliced; $k++){
		my $tr1=$list1Spliced[$k];
		print STDERR "DEBUG: tr1 $transcript_index{$tr1}\n" if $debug;
		if(exists ($transcript_to_container{$tr1})){
			print STDERR "DEBUG: tr1 IS CONTAINED in $transcript_index{$transcript_to_container{$tr1}}\n" if $debug;
			next;
		}
    @{$container_original_ends{$tr1}}=(${$transcript_exons{$tr1}}[0][1], ${$transcript_exons{$tr1}}[-1][2]);
		next unless ($supportCount{$tr1} >= $minReadSupport);

		NEXTTR2: for (my $l=$k+1; $l<=$#list1Spliced; $l++){
			my $tr2=$list1Spliced[$l];
			print STDERR "DEBUG:  tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2}\n" if $debug;
			next unless ($supportCount{$tr2} >= $minReadSupport);
			if(exists ($transcript_to_container{$tr2})){
				next;
			}

			if(${$transcript_exons{$tr1}}[-1][2] - ${$transcript_exons{$tr2}}[0][1]>0){
				if( $transcript_strand{$tr1} == $transcript_strand{$tr2}){
					#are intron chains compatible? (ie is tr2 a subset or equal to tr1?
					my $intronTr1Index=0;
					my $countIntronsTr2MatchedToTr1=-1;
					my $firstTr1IntronMatchedToTr2;
					my $lastTr1IntronMatchedToTr2;
					for (my $j=0; $j <= $#{$transcript_introns{$tr2}}; $j++){
						print STDERR "DEBUG:   intronTr2 #$j: ${${$transcript_introns{$tr2}}[$j]}[1], ${${$transcript_introns{$tr2}}[$j]}[2]\n" if $debug;
						for (my $i=$intronTr1Index; $i<= $#{$transcript_introns{$tr1}}; $i++){
							print STDERR "DEBUG:    intronTr1 #$i: ${${$transcript_introns{$tr1}}[$i]}[1], ${${$transcript_introns{$tr1}}[$i]}[2]\n" if $debug;


							### BEGIN this should be a subroutine but repeated subroutine calls are too expensive in perl
							### my $overlap=overlap(${$intronTr1}[1], ${$intronTr1}[2], ${$intronTr2}[1], ${$intronTr2}[2]);
							#my $start1=${${$transcript_introns{$tr1}}[$i]}[1];
							#my $stop1=${${$transcript_introns{$tr1}}[$i]}[2];
							#my $start2=${${$transcript_introns{$tr2}}[$j]}[1];
							#my $stop2=${${$transcript_introns{$tr2}}[$j]}[2];
							my $overlap;
							my $start2minusstop1=${${$transcript_introns{$tr2}}[$j]}[1]-${${$transcript_introns{$tr1}}[$i]}[2];
							my $start2minusstart1=${${$transcript_introns{$tr2}}[$j]}[1]-${${$transcript_introns{$tr1}}[$i]}[1];
							my $stop2minusstart1=${${$transcript_introns{$tr2}}[$j]}[2]-${${$transcript_introns{$tr1}}[$i]}[1];
							my $stop2minusstop1=${${$transcript_introns{$tr2}}[$j]}[2]-${${$transcript_introns{$tr1}}[$i]}[2];
							if( ( $stop2minusstart1>=0  && $stop2minusstop1 <=0 ) || ($start2minusstart1 >=0 && $start2minusstop1 <=0) || ($start2minusstart1 <= 0 && $stop2minusstop1 >= 0)){
								$overlap=1;
							}
							else{
								if($stop2minusstart1<0){
									$overlap=-1
								}
								elsif($start2minusstop1>0){
									$overlap=-2
								}
							}
							### # 1  : overlap
		                        # -1 : 2 upstream of 1
		                        # -2 : 2 downstream of 1
							### END this should be a subroutine but repeated subroutine calls are too expensive in perl



							print STDERR "DEBUG:     overlap $overlap\n" if $debug;
							if($overlap == 1){
								if(${${$transcript_introns{$tr2}}[$j]}[1] == ${${$transcript_introns{$tr1}}[$i]}[1] && ${${$transcript_introns{$tr2}}[$j]}[2] == ${${$transcript_introns{$tr1}}[$i]}[2]){
									$countIntronsTr2MatchedToTr1++;
									if($j==0){
										$firstTr1IntronMatchedToTr2=$i;
									}
									if($j==$#{$transcript_introns{$tr2}}){
										$lastTr1IntronMatchedToTr2=$i;
									}
									$intronTr1Index=$i+1; #skip directly to next $tr1 intron at the next round (next $tr2 intron)
									if($intronTr1Index > $#{$transcript_introns{$tr1}} #we've reached the last intron of tr1.
										&& $j < $#{$transcript_introns{$tr2}} #we've not reached the last intron of tr2.
										&& ${$transcript_introns{$tr1}}[0][1] == ${$transcript_introns{$tr2}}[0][1] #tr1 and tr2's respective intron chains start at the same coord
										&& ${$transcript_introns{$tr2}}[$j+1][1] > ${$transcript_exons{$tr1}}[-1][2] - $exonOverhangTolerance){ #tr1's last exon does not overhang too much inside tr2's next intron
										#Transfer remaining exons of tr2 to tr1, if they're compatible
										$lastTr1IntronMatchedToTr2=$i;
										print STDERR "DEBUG:     reached last tr1 intron\n" if $debug;
										if(checkIntronExonOverlap($tr1,$tr2,$firstTr1IntronMatchedToTr2,$lastTr1IntronMatchedToTr2) ==0){

											$container_to_transcripts{$tr1}{$tr2}=undef ;
											$transcript_to_container{$tr2}=$tr1;
											transferRightExonsIntrons($tr1,$tr2,$j);
											adjustContainerEnds($tr1,$tr2);
											if (exists $container_to_transcripts{$tr2}){
												foreach my $tr3 (keys %{$container_to_transcripts{$tr2}}){
													$transcript_to_container{$tr3}=$tr1;
													$container_to_transcripts{$tr1}{$tr3}=undef ;
												}
												delete($container_to_transcripts{$tr2});
											}

											print STDERR "DEBUG:     tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2} MATCH 2\n" if $debug;
											next NEXTTR2

										}
										else{
											print STDERR "DEBUG:   tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2} INCOMP 2\n" if $debug;
											next NEXTTR2;
										}
									}
							 		last;
								}
								else{ # introns overlap but don't exactly match, give up current tr2
									print STDERR "DEBUG:    tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2} INCOMP 3\n" if $debug;
									next NEXTTR2;
								}
							}
							elsif($overlap == -1){ #intron1 is downstream of intron2
								if($countIntronsTr2MatchedToTr1 == $#{$transcript_introns{$tr2}}){ # all tr2 introns have found a match in tr1
									print STDERR "DEBUG:    YES 1" if $debug;
									if(checkIntronExonOverlap($tr1,$tr2,$firstTr1IntronMatchedToTr2,$lastTr1IntronMatchedToTr2) ==0){
										$container_to_transcripts{$tr1}{$tr2}=undef ;
										$transcript_to_container{$tr2}=$tr1;
										adjustContainerEnds($tr1,$tr2);
										if (exists $container_to_transcripts{$tr2}){
											foreach my $tr3 (keys %{$container_to_transcripts{$tr2}}){
												$transcript_to_container{$tr3}=$tr1;
												$container_to_transcripts{$tr1}{$tr3}=undef ;
											}
											delete($container_to_transcripts{$tr2});
										}

										print STDERR "DEBUG:   tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2} MATCH 3\n" if $debug;
										next NEXTTR2;
									}
									else{
										print STDERR "DEBUG:   tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2} INCOMP 4\n" if $debug;
										next NEXTTR2;
									}
								}
								else{

									if($i==0 && $j==0){ #first intron of tr1 is downstream of first intron of tr2. Since transcripts are sorted by position of first intron, we can skip the rest of the tr2 list
										print STDERR "DEBUG:   tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2} INCOMP 5\n" if $debug;
										next NEXTTR1;
									}
									else{
										print STDERR "DEBUG:   tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2} INCOMP 6\n" if $debug;
										next NEXTTR2;
									}
								}
							}
							elsif($overlap == -2){ #intron1 is upstream of intron2
								if($i == $#{$transcript_introns{$tr1}} || $j > 0 ){ #last intron of tr1 is upstream of first intron of tr2, i.e. tr1 and tr2 are incompatible
									print STDERR "DEBUG:   tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2} INCOMP 7\n" if $debug;
									next NEXTTR2;
								}

							}
						}
					}
					if($countIntronsTr2MatchedToTr1 == $#{$transcript_introns{$tr2}}){ #fully identical or contained intron chain
						print STDERR "DEBUG:    YES 2\n" if $debug;
						if(checkIntronExonOverlap($tr1,$tr2,$firstTr1IntronMatchedToTr2,$lastTr1IntronMatchedToTr2) ==0){
							$container_to_transcripts{$tr1}{$tr2}=undef ;
							$transcript_to_container{$tr2}=$tr1;
							adjustContainerEnds($tr1,$tr2);
							if (exists $container_to_transcripts{$tr2}){
								foreach my $tr3 (keys %{$container_to_transcripts{$tr2}}){
									$transcript_to_container{$tr3}=$tr1;
									$container_to_transcripts{$tr1}{$tr3}=undef ;
								}
								delete($container_to_transcripts{$tr2});
							}
							print STDERR "DEBUG:    tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2} MATCH 4\n" if $debug;
							next;
						}
					}
					else{
						print STDERR "DEBUG:   tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2} INCOMP 9\n" if $debug;

					}

				}
			}

			else{ #tr2 is downstream of tr1. Since trs are sorted by start position, we can safely skip to next tr1
			 	print STDERR "DEBUG:   tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2} INCOMP 10\n" if $debug;
			 	last;
			}

		}
	}
	%transcript_introns=();

	@list1Mono=sort { ${$transcript_exons{$a}}[0][1] <=> ${$transcript_exons{$b}}[0][1] or ${$transcript_exons{$a}}[0][2] <=> ${$transcript_exons{$b}}[0][2] } @list1Mono;
	# compute read support for each input monoexonic transcript/read:
	print STDERR "Calculating monoexonic read support...\n" if $debug;
	NEXTTR1: for (my $k=0; $k<=$#list1Mono; $k++){
		my $tr1=$list1Mono[$k];
		print STDERR "DEBUG: tr1 $transcript_index{$tr1}\n" if $debug;
		if(exists ($transcript_to_container{$tr1})){
			print STDERR "DEBUG: tr1 IS CONTAINED in $transcript_index{$transcript_to_container{$tr1}}\n" if $debug;
			next;
		}
    @{$container_original_ends{$tr1}}=(${$transcript_exons{$tr1}}[0][1], ${$transcript_exons{$tr1}}[-1][2]);
		NEXTTR2: for (my $l=$k+1; $l<=$#list1Mono; $l++){
			my $tr2=$list1Mono[$l];

			print STDERR "DEBUG:  tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2}\n" if $debug;
			if(exists ($transcript_to_container{$tr2})){
				next;
			}
			if(${$transcript_exons{$tr1}}[0][2] - ${$transcript_exons{$tr2}}[0][1]>=0){
				if( ($transcript_strand{$tr1} == $transcript_strand{$tr2})
				   &&
				   (${$transcript_exons{$tr2}}[0][1] >= ${$transcript_exons{$tr1}}[0][1] - $transcriptEndFuzziness
				      &&
				      ${$transcript_exons{$tr2}}[0][1] <= ${$transcript_exons{$tr1}}[0][1] + $transcriptEndFuzziness)
				      &&
				      (${$transcript_exons{$tr2}}[-1][2] >= ${$transcript_exons{$tr1}}[-1][2] - $transcriptEndFuzziness
				      &&
				      ${$transcript_exons{$tr2}}[-1][2] <= ${$transcript_exons{$tr1}}[-1][2] + $transcriptEndFuzziness) ) {

					$container_to_transcripts{$tr1}{$tr2}=undef ;
					$transcript_to_container{$tr2}=$tr1;
					adjustContainerEnds($tr1,$tr2);
					if (exists $container_to_transcripts{$tr2}){
						foreach my $tr3 (keys %{$container_to_transcripts{$tr2}}){
							$transcript_to_container{$tr3}=$tr1;
							$container_to_transcripts{$tr1}{$tr3}=undef ;
						}
						delete($container_to_transcripts{$tr2});
					}

					$supportCount{$tr1}++;
					$container_to_supporting_transcripts{$tr1}{$tr2}=1;
					print STDERR "DEBUG:     tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2} FULLMATCH 2\n" if $debug;
				}
			}
			else{ #tr2 is downstream of tr1. Since trs are sorted by start position, we can safely skip to next tr1
		 		print STDERR "DEBUG:   tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2} INCOMP 10\n" if $debug;
		 		last;
			}

		}

	}

	print STDERR "Merging monoexonic reads...\n" if $debug;
	NEXTTR1: for (my $k=0; $k<=$#list1Mono; $k++){
		my $tr1=$list1Mono[$k];
		print STDERR "DEBUG: tr1 $transcript_index{$tr1}\n" if $debug;
		if(exists ($transcript_to_container{$tr1})){
			print STDERR "DEBUG: tr1 IS CONTAINED in $transcript_index{$transcript_to_container{$tr1}}\n" if $debug;
			next;
		}
    @{$container_original_ends{$tr1}}=(${$transcript_exons{$tr1}}[0][1], ${$transcript_exons{$tr1}}[-1][2]);
		print STDERR "DEBUG: tr1 read support: $supportCount{$tr1}\n" if $debug;
		next unless ($supportCount{$tr1} >= $minReadSupport);
		NEXTTR2: for (my $l=$k+1; $l<=$#list1Mono; $l++){
			my $tr2=$list1Mono[$l];
			print STDERR "DEBUG:  tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2}\n" if $debug;
			print STDERR "DEBUG:   tr2 read support: $supportCount{$tr2}\n" if $debug;
			next unless ($supportCount{$tr2} >= $minReadSupport);
			if(exists ($transcript_to_container{$tr2})){
				print STDERR "DEBUG:   tr2 is contained in $transcript_to_container{$tr2}\n" if $debug;
				next;
			}
			if(${$transcript_exons{$tr1}}[0][2] - ${$transcript_exons{$tr2}}[0][1]>=0){
				if( $transcript_strand{$tr1} == $transcript_strand{$tr2}){
					$container_to_transcripts{$tr1}{$tr2}=undef ;
					$transcript_to_container{$tr2}=$tr1;
					adjustContainerEnds($tr1,$tr2);
					if (exists $container_to_transcripts{$tr2}){
						foreach my $tr3 (keys %{$container_to_transcripts{$tr2}}){
							$transcript_to_container{$tr3}=$tr1;
							$container_to_transcripts{$tr1}{$tr3}=undef ;
						}
						delete($container_to_transcripts{$tr2});
					}

					print STDERR "DEBUG:  tr1 $transcript_index{$tr1} overlaps tr2 $transcript_index{$tr2}\n" if $debug;
				}

			}
			else{
				print STDERR "DEBUG:   tr1 $transcript_index{$tr1} vs tr2 $transcript_index{$tr2} INCOMP 11\n" if $debug;
				last;

			}
		}

	}


  ################################
  ################################
	##     print all to GTF       ##
  ################################
  ################################

	foreach my $tr (@list1Spliced){
		next unless ($supportCount{$tr} >= $minReadSupport);
		if (exists ($container_to_transcripts{$tr})){ # i.e. if tr is a container
			my @trList= sort { ${$transcript_exons{$b}}[-1][2] - ${$transcript_exons{$b}}[0][1] <=> ${$transcript_exons{$a}}[-1][2] - ${$transcript_exons{$a}}[0][1] } (keys (%{$container_to_transcripts{$tr}}), $tr); # sort by genomic length
			my @longest=();
			my $longestTMlength=0;
			for (my $i=0;$i<=$#trList;$i++){
				my $t=$trList[$i];
				if($i==0){
					$longestTMlength=(${$transcript_exons{$t}}[-1][2] - ${$transcript_exons{$t}}[0][1]);
				}
				if(${$transcript_exons{$t}}[-1][2] - ${$transcript_exons{$t}}[0][1] == $longestTMlength){
					push(@longest, $t);
				}
			}
      my $trMatureLength=0;
      for (my $i=0; $i<=$#{$transcript_exons{$tr}}; $i++){
        $trMatureLength+= (${${$transcript_exons{$tr}}[$i]}[2] - ${${$transcript_exons{$tr}}[$i]}[1])+1;
      }
			#my @transcriptsList=();
			#foreach my $t (@trList){
			#	push(@transcriptsList, $t)
			#}
			my %uniqsupporting_longest=();
			my @supporting_longest=();
			#print STDERR "longest:".join(",", @longest)."\n";
			foreach my $longest (@longest){
				foreach my $supporting_longest (keys %{$container_to_supporting_transcripts{$longest}}, $longest){
					#print STDERR "supporting: $supporting\n";
					$uniqsupporting_longest{$supporting_longest}=1;
				}
			}
			#print STDERR Dumper \%uniqSupporting;

			foreach my $t (keys %uniqsupporting_longest){
				#print STDERR "$t\n";
				push(@supporting_longest, $t);
			}
			if ($debug){
				print STDERR "print records for ".$transcript_index{$tr}."\n";
				print STDERR "\ttrList:\n";
				foreach my $id (@trList){
					print STDERR "\t ".$transcript_index{$id}."\n";
				}
			}

      #compute distance of all contained reads' TSSs/TTSs to TM's TSS/TTS
      my ($refDistsToLeft, $refDistsToRight) = endDistances($tr, \@trList);
      #print STDERR "Dists to left: ".join(",", @{$refDistsToLeft})."\n";
      #print STDERR "Dists to right: ".join(",", @{$refDistsToRight})."\n";

			printGTF($tr, 0, \@trList, \@longest, \@supporting_longest, 1, $refDistsToLeft, $refDistsToRight, $trMatureLength);
		}
		else{ #tr is not a container
			unless(exists ($transcript_to_container{$tr})) { # i.e. if tr is not a container, and is not contained
			my @trList=($tr);
      my @tmp1=('0');
      my $trMatureLength=0;
      for (my $i=0; $i<=$#{$transcript_exons{$tr}}; $i++){
        $trMatureLength+= (${${$transcript_exons{$tr}}[$i]}[2] - ${${$transcript_exons{$tr}}[$i]}[1])+1;
      }

      #print STDERR "Dists to left: ".join(",", @tmp1)."\n";
      #print STDERR "Dists to right: ".join(",", @tmp1)."\n";
			if ($debug){
				print STDERR "print records for ".$transcript_index{$tr}."\n";
				print STDERR "\ttrList:\n";
				foreach my $id (@trList){
					print STDERR "\t ".$transcript_index{$id}."\n";
				}
			}
				printGTF($tr, 0, \@trList, \@trList, \@trList, 1, \@tmp1, \@tmp1, $trMatureLength);
			}
		}
	}


	foreach my $tr (@list1Mono){
		next unless ($supportCount{$tr} >= $minReadSupport);
		if (exists ($container_to_transcripts{$tr})){ # i.e. if tr is a container
			my @trList= sort { ${$transcript_exons{$b}}[-1][2] - ${$transcript_exons{$b}}[0][1] <=> ${$transcript_exons{$a}}[-1][2] - ${$transcript_exons{$a}}[0][1] } (keys (%{$container_to_transcripts{$tr}}), $tr); # sort by genomic length
			if ($#trList>=$minReadSupport-1){
				my @longest=();
				my $longestTMlength=0;
				for (my $i=0;$i<=$#trList;$i++){
					my $t=$trList[$i];
					if($i==0){
						$longestTMlength=${$transcript_exons{$t}}[-1][2] - ${$transcript_exons{$t}}[0][1];
					}
					if(${$transcript_exons{$t}}[-1][2] - ${$transcript_exons{$t}}[0][1] == $longestTMlength){
						push(@longest, $t);
					}
				}
        my $trMatureLength=0;
        for (my $i=0; $i<=$#{$transcript_exons{$tr}}; $i++){
          $trMatureLength+= (${${$transcript_exons{$tr}}[$i]}[2] - ${${$transcript_exons{$tr}}[$i]}[1])+1;
        }

				#my @transcriptsList=();
				#foreach my $t (@trList){
				#	push(@transcriptsList, $t)
				#}
				my %uniqsupporting_longest=();
				my @supporting_longest=();
				foreach my $longest (@longest){
					foreach my $supporting_longest (keys %{$container_to_supporting_transcripts{$longest}}, $longest){
						$uniqsupporting_longest{$supporting_longest}=1;
					}
				}
				foreach my $t (keys %uniqsupporting_longest){
					push(@supporting_longest, $t);
				}
      #compute distance of all contained reads' TSSs/TTSs to TM's TSS/TTS
      my ($refDistsToLeft, $refDistsToRight) = endDistances($tr, \@trList);
      #print STDERR "Dists to left: ".join(",", @{$refDistsToLeft})."\n";
      #print STDERR "Dists to right: ".join(",", @{$refDistsToRight})."\n";
				printGTF($tr, 0, \@trList, \@longest, \@supporting_longest, 0, $refDistsToLeft, $refDistsToRight, $trMatureLength);
			}

		}
		else{ #tr is not a container
			unless(exists ($transcript_to_container{$tr})) { # i.e. if tr is not a container, and is not contained
				my @trList=($tr);
        my @tmp1=('0');
        my $trMatureLength=0;
        for (my $i=0; $i<=$#{$transcript_exons{$tr}}; $i++){
          $trMatureLength+= (${${$transcript_exons{$tr}}[$i]}[2] - ${${$transcript_exons{$tr}}[$i]}[1])+1;
        }

        #print STDERR "Dists to left: ".join(",", @tmp1)."\n";
        #print STDERR "Dists to right: ".join(",", @tmp1)."\n";
				printGTF($tr, 0, \@trList, \@trList, \@trList, 0, \@tmp1, \@tmp1, $trMatureLength);

			}
		}
	}

}

print STDERR "Done.\n";



sub buildContig{
	my $trA=$_[0];
	my $feature_to_feature=$_[1];
	my $feature_to_contig=$_[2];
	my $contigNumber=$_[3];

	unless(exists ${$feature_to_contig}{$trA}){
		${$feature_to_contig}{$trA}=$contigNumber;
	}
	foreach my $trB (keys %{${$feature_to_feature}{$trA}}){
		unless( exists ${$feature_to_contig}{$trB} ){
			buildContig($trB, $feature_to_feature, $feature_to_contig, $contigNumber);
		}
	}
}

sub printGTF{
	my $transcript_id=$_[0];
	my $score=$_[1];
	my @contains=@{$_[2]};
	my @longest=@{$_[3]};
	my @supporting_longest=@{$_[4]};
	my $supporting_longest_count= scalar @supporting_longest;
	my $contains_count= scalar @contains;
	my $rpm=$contains_count/$million_read_count;
	my $flrpm=$supporting_longest_count/$million_read_count;
	my $spliced_bool=$_[5];
  my @distsToLeft=@{$_[6]};
  my @distsToRight=@{$_[7]};
  my $length=$_[8];
  my @distsToFivePend=();
  my @distsToThreePend=();
	my @contains_id=();
	foreach my $id (@contains){
		push(@contains_id, $transcript_index{$id});
	}

	my @longest_id=();
	foreach my $id (@longest){
		push(@longest_id, $transcript_index{$id});
	}

	my @supporting_longest_id=();
	foreach my $id (@supporting_longest){
		push(@supporting_longest_id, $transcript_index{$id});
	}
	my $contains=join(",", @contains_id);
	my $longest=join(",", @longest_id);
	my $supporting_longest=join(",", @supporting_longest_id);
  my $strand='';
  if(${$transcript_exons{$transcript_id}}[0][3] == -1){
    $strand='-';
    @distsToFivePend=@distsToRight;
    @distsToThreePend=@distsToLeft;
  }
  elsif (${$transcript_exons{$transcript_id}}[0][3] == 1){
    $strand='+';
    @distsToFivePend=@distsToLeft;
    @distsToThreePend=@distsToRight;
  }
  elsif (${$transcript_exons{$transcript_id}}[0][3] == 0){
    $strand='.';
    @distsToFivePend=@distsToLeft;
    @distsToThreePend=@distsToRight;
  }
  else{
    die;
  }
  my $distsToFivePendString=join(",", @distsToFivePend);
  my $distsToThreePendString=join(",", @distsToThreePend);
  my @metaDistsToFivePend=();
  my @metaDistsToThreePend=();
  foreach my $dist (@distsToFivePend){
    push(@metaDistsToFivePend, $dist/$length);
  }
  foreach my $dist (@distsToThreePend){
    push(@metaDistsToThreePend, 1 - ($dist/$length));
  }
  my $metaDistsToFivePendString=join(",", @metaDistsToFivePend);
  my $metaDistsToThreePendString=join(",", @metaDistsToThreePend);


	foreach my $exon (@{$transcript_exons{$transcript_id}}){
		my $tmId=makeTmId($transcript_id);
		print "${$exon}[0]\ttmerge\texon\t${$exon}[1]\t${$exon}[2]\t$score\t$strand\t.\tgene_id \"$tmId\"; transcript_id \"$tmId\"; contains \"$contains\"; contains_count \"$contains_count\"; 3p_dists_to_3p \"$distsToThreePendString\"; 5p_dists_to_5p \"$distsToFivePendString\"; flrpm \"$flrpm\"; longest \"$longest\"; longest_FL_supporters \"$supporting_longest\"; longest_FL_supporters_count \"$supporting_longest_count\"; mature_RNA_length \"$length\"; meta_3p_dists_to_5p \"$metaDistsToThreePendString\"; meta_5p_dists_to_5p \"$metaDistsToFivePendString\"; rpm \"$rpm\"; spliced \"$spliced_bool\";\n";
	}
}

sub endDistances{
  my $container=$_[0];
  my @contains=@{$_[1]};
  my @distsToLeft=();
  my @distsToRight=();
  #print STDERR $transcript_index{$container}."\n";
  foreach my $trContained (@contains){
    #print STDERR "\t$transcript_index{$trContained}\n";
    my $containedLeftStart;
    my $containedRightEnd;
    if($trContained == $container){
      $containedLeftStart=${$container_original_ends{$container}}[0];
      $containedRightEnd=${$container_original_ends{$container}}[1];
    }
    else{
      $containedLeftStart=${$transcript_exons{$trContained}}[0][1];
      $containedRightEnd=${$transcript_exons{$trContained}}[-1][2];
    }
    #print STDERR "\t$containedLeftStart - $containedRightEnd\n";

    my $distToLeft=0;
    my $intronLeftSubstract=0;
    my $distToRight=0;
    my $intronRightSubstract=0;

    #compute distance to left end of container
    for (my $i=0; $i<=$#{$transcript_exons{$container}}; $i++){
      if($i>0){
        $intronLeftSubstract+=${${$transcript_exons{$container}}[$i]}[1] - ${${$transcript_exons{$container}}[$i-1]}[2];
        $distToLeft=(${${$transcript_exons{$container}}[$i]}[1] - ${${$transcript_exons{$container}}[0]}[1]) - $intronLeftSubstract;
      }
      #print STDERR "\t i: $i ${${$transcript_exons{$container}}[$i]}[1] ${${$transcript_exons{$container}}[$i]}[2]\n";
      #print STDERR "\t i: $i distToLeft: $distToLeft\n";
      if($containedLeftStart >= (${${$transcript_exons{$container}}[$i]}[1] -1) - $exonOverhangTolerance && $containedLeftStart <= ${${$transcript_exons{$container}}[$i]}[2]){
        my $lastDist=$containedLeftStart - ${${$transcript_exons{$container}}[$i]}[1] ;
        $lastDist=0 if $lastDist<0; #account for exonOverhangTolerance
        $distToLeft+=$lastDist ;
        last;
      }
      elsif ($containedLeftStart < ${${$transcript_exons{$container}}[$i]}[1]){
        die "ERROR: Program died due to a bug (in endDistances subroutine: $containedLeftStart < ${${$transcript_exons{$container}}[$i]}[1])\nsorry. Please contact author.\n";
      }
    }

    #compute distance to right end of container
    for (my $i=$#{$transcript_exons{$container}}; $i>=0; $i--){
      if($i<$#{$transcript_exons{$container}}){
        $intronRightSubstract+=${${$transcript_exons{$container}}[$i+1]}[1] - ${${$transcript_exons{$container}}[$i]}[2];
        $distToRight=(${${$transcript_exons{$container}}[$#{$transcript_exons{$container}}]}[2] - ${${$transcript_exons{$container}}[$i]}[2]) - $intronRightSubstract;
        #$distToLeft=$distToLeft + (${${$transcript_exons{$container}}[$i]}[1] - ${${$transcript_exons{$container}}[$i-1]}[1]);
        #$distToLeft= $distToLeft - (${${$transcript_exons{$container}}[$i]}[1] - ${${$transcript_exons{$container}}[$i-1]}[2]);
      }
      #print STDERR "\t i: $i distToRight: $distToRight\n";
      if($containedRightEnd <= (${${$transcript_exons{$container}}[$i]}[2] +1) + $exonOverhangTolerance && $containedRightEnd >= ${${$transcript_exons{$container}}[$i]}[1]){
        my $lastDist=${${$transcript_exons{$container}}[$i]}[2] - $containedRightEnd ;
        $lastDist=0 if $lastDist<0; #account for exonOverhangTolerance
        $distToRight+=$lastDist ;
        last;
      }
      elsif ($containedRightEnd > ${${$transcript_exons{$container}}[$i]}[2]){
        die "ERROR: Program died due to a bug (in endDistances subroutine: $containedRightEnd > ${${$transcript_exons{$container}}[$i]}[2])\n, sorry. Please contact author.\n";
      }
    }
    #print STDERR "\tFinal distToLeft: $distToLeft\n";
    #print STDERR "\tFinal distToRight: $distToRight\n";
    push(@distsToLeft, $distToLeft);
    push(@distsToRight, $distToRight);

  }
  return (\@distsToLeft, \@distsToRight);
}



sub checkIntronExonOverlap{
#verify that when two A and B transcripts have compatible intron chains, there is no terminal exon/intron overlap
#returns 1 if any overlap found, 0 otherwise
	my $trA=$_[0];
	my $trB=$_[1];
	my $firstTrAIntronMatchedToTrB=$_[2];
	my $lastTrAIntronMatchedToTrB=$_[3];
	unless (defined $firstTrAIntronMatchedToTrB){
		die "ERROR: firstTrAIntronMatchedToTrB is undefined for $transcript_index{$trA} / $transcript_index{$trB}\n";
	}
	unless (defined $lastTrAIntronMatchedToTrB){
		die "ERROR: lastTrAIntronMatchedToTrB is undefined for $transcript_index{$trA} / $transcript_index{$trB}\n";

	}
		my @exons2=(${$transcript_exons{$trB}}[0],${$transcript_exons{$trB}}[-1]); #only terminal exons
		for (my $i=0; $i<=$#exons2;$i++) {
			my $exon2=$exons2[$i];
			if($i==0){ #first exon of trB
				unless ($firstTrAIntronMatchedToTrB ==0) { # if the first intron of trA matched to trB is trA[0], no need to look for upstream introns on trA

					my $intron1=${$transcript_introns{$trA}}[$firstTrAIntronMatchedToTrB-1];
					#my $overlap=overlap(${$intron1}[1], ${$intron1}[2], ${$exon2}[1], ${$exon2}[2]);
					if( ${$exon2}[1] <= ${$intron1}[2] - $exonOverhangTolerance){
						return 1;
					}
				}
			}
			elsif($i==1){ #last exon of trB
				unless ($lastTrAIntronMatchedToTrB == $#{$transcript_introns{$trA}} ){ # if the last intron of trA matched to trB is trA[-1], no need to look for downstream introns on trA
					my $intron1=${$transcript_introns{$trA}}[$lastTrAIntronMatchedToTrB+1];
					#my $overlap=overlap(${$intron1}[1], ${$intron1}[2], ${$exon2}[1], ${$exon2}[2]);
					if(${$exon2}[2] >= ${$intron1}[1] + $exonOverhangTolerance){
						return 1;
					}
				}
			}
			else{
				die;
			}

		}

	return 0;
}



sub makeTmId{
	my $id=$_[0];
	my @newId=split("", $id);
	my @prepend=(join("",split("", $tmPrefix)),'T','M','_');
	my $totalLength=($#prepend+13)-length($id);
	for (my $i=$#prepend+1;$i<$totalLength; $i++){
		$prepend[$i]=0;
	}
	unshift(@newId, @prepend);
	return(join("",@newId))
}


sub adjustContainerEnds{
	my $tr1=$_[0];
	my $tr2=$_[1];
	if(${$transcript_exons{$tr2}}[0][1] < ${$transcript_exons{$tr1}}[0][1]){
		${$transcript_exons{$tr1}}[0][1] = ${$transcript_exons{$tr1}}[0][1]
	}
	if(${$transcript_exons{$tr2}}[-1][2] > ${$transcript_exons{$tr1}}[-1][2]){
		${$transcript_exons{$tr1}}[-1][2] = ${$transcript_exons{$tr2}}[-1][2]
	}
}

sub transferRightExonsIntrons{
	my $tr1=$_[0]; #container
	my $tr2=$_[1]; #content
	my $tr2IntronIndex=$_[2];
	my $leftmostTr2ExonToTransfer=$tr2IntronIndex+1;
	pop @{$transcript_exons{$tr1}};
	for (my $m=$leftmostTr2ExonToTransfer; $m <= $#{$transcript_exons{$tr2}} ;$m++) {
		push (@{$transcript_exons{$tr1}}, \@{${$transcript_exons{$tr2}}[$m]} );
	}
	for (my $m=$tr2IntronIndex+1; $m <= $#{$transcript_introns{$tr2}} ;$m++) {
		push (@{$transcript_introns{$tr1}}, \@{${$transcript_introns{$tr2}}[$m]} );
	}

}
