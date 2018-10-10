#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

my $program = 'unique';					#name of script

my $seq_recog = "acgtn\?:~-";					# MOST AMBIGUITY CODES NOT RECOGNISED FOR MOST SEQUENCES

my @changes = qw(A->T A->C A->G C->A C->G C->T G->A G->T G->C T->A T->C T->G);
my @bases = qw(A C G T);

my %parameters;							#input and output file names
my ($datafile, $prefix,$outfile,$fastaindex,$lengthfile,$nonuniquefile,$insertionfile,$ambancfile,$hapfile,$hapseqfile);		
my (@names, %seq);
my $alignment_length = 0;
my $refpattern = 'ORF';
my $insthreshold = 0;
my $getuniqueatshared = "no";
my $ucscpos = "no";
my $includeref = "no";
my $hotspotlength = 7;

my (@refseqnames, @neighbour_changes, @neighbour_bases);


($prefix = $program) =~ s/(.*?).pl/$1/; 					# default output prefix

getopts('i:n:m:o:r:I:H:fhtsCNcuRMk',\%parameters);


if (exists($parameters{"i"})) { $datafile = $parameters{"i"}; }
if (exists($parameters{"n"})) { $datafile = $parameters{"n"}; }
if (exists($parameters{"m"})) { $datafile = $parameters{"m"}; }
if (exists($parameters{"o"})) { $prefix = $parameters{"o"}; }
if (exists($parameters{"r"})) { $refpattern = $parameters{"r"}; }
if (exists($parameters{"R"})) { $includeref = "yes"; }
if (exists($parameters{"I"})) { 
	if (defined $parameters{"I"}) { $insthreshold = $parameters{"I"}; }
}										# OPTIONAL:
if (exists($parameters{"u"})) { $ucscpos = "yes"; }				#default is to output as 1..100
if (exists $parameters{"f"}) { $fastaindex = $prefix."_withindex.fa"; }		#default fastaindex filename
if (exists $parameters{"s"}) { $nonuniquefile = $prefix."_nonunique.fa"; } 	#default nonunique filename
if (exists $parameters{"N"}) { $getuniqueatshared = "yes"; } 			#count unique point subs at shared sites (manual labour required by user for final calls)
if (exists $parameters{"H"}) { $hotspotlength = $parameters{"H"}; } 	
if (exists $parameters{'c'}) {
	foreach my $change (@changes) { foreach (0,1,2,"?") { push (@neighbour_changes, "$_$change"); } }
	unless (exists $parameters{"N"}) { foreach (0,1,2,"?") { push (@neighbour_changes, "$_?"); } }

	foreach my $base (@bases) { foreach (0,1,2,"?") { push (@neighbour_bases, "$_$base"); } }	
}

										# NOT OPTIONAL:
$outfile = $prefix.".subs";							#default outfile name
$lengthfile = $prefix.".lengths";						#default length filename
$insertionfile = $prefix.".ins";				 		#default insertionfile name
$ambancfile = $prefix."ambanc.xls";				 		#default ambiguous ancestors filename
$hapfile =  $prefix.".haps";							#default haplotype filename	
$hapseqfile =  $prefix.".haps.gfa";							#default haplotype alignment filename

print "\n\t< $program > - Counting unique subs in an alignment\n\n";

unless ((exists $parameters{"i"}) || (exists $parameters{"n"}) || (exists $parameters{"m"})) {
	print "\n USAGE: $program -i <datafile>\n";
	print   " OR:    $program -n <datafile>\n"; 
	print   " OR:    $program -m <datafile>\n\n"; 
	print   "    -i\tdatafile: sequence alignment in fasta format\n";
	print 	"    -n\tdatafile: sequence alignment in nexus format\n";
	print 	"    -m\tdatafile: sequence alignment in mase format\n";
	print   "    -o\tprefix for output files [default: $prefix]\n";	
	print   "    -r\tpattern for reference sequence names [default: 1st sequence]\n";
	print   "      \tto define open reading frame(s).\n"; 
	print   "    -k\tcreate a majority consensus to use as reference sequence\n";
	print   "      \tto define open reading frame(s) or direction of indels.\n"; 
	print   "    -R\tinclude reference sequence(s) when deciding unique sites\n";
	print   "      \tunique changes in reference sequences are not counted.\n"; 
	print   "    -M\tmerge identical sequences into a single haplotype\n";
	print   "    -t\tdo not count 5' truncation as deletions [e.g. for nonLTRs].\n"; 
	print   "    -C\ttreat complex events as separate insertions and deletions.\n"; 
	print   "    -N\tcount unique point subs at non-unique sites.\n"; 	
	print   "    -H\tminimum length of homopolymer runs in reference seq to report [$hotspotlength bp]\n"; 
	print   "    -c\tshow nearest neighbour context for unique point subs.\n";
	print   "    -I\tminimum length of Insertions to print to [prefix]_ins.fa\n"; 
	print   "    -s\tprint sites with non-unique point subs to fasta file [prefix]_nonunique.fa.\n"; 
	print   "    -f\tprint fasta alignment with index as 1st seq.\n"; 
	print   "    -u\treport nucleotide positions in ucsc format ie 0 based start, 1 based end [1..100].\n"; 
	print   "    -h\thide headers for output tables.\n\n"; 

	print   "NOTE: this script does not recognize ambiguity codes e.g. RYWS etc. Convert these to Ns before using. fastaLC2n.pl with the -A option will do this. Also this script needs too much memory for genome alignments. Use al2variable.pl to reduce the size of your alignment to variable sites only in that case.\n\n"; 

	exit;
}

print "\ncommand: $program ";
foreach my $flag (sort keys %parameters) { 
	print "-$flag "; 
	if (defined $parameters{$flag}) { print "$parameters{$flag} "; }
}
print "\n\n";

## READ IN SEQUENCES FROM FILE 
#
my $seq_ref;

if (exists $parameters{'n'}) { ($seq_ref,$alignment_length,@names) = nexus($datafile); }	# NEXUS FORMAT
elsif (exists $parameters{'m'}) { ($seq_ref,$alignment_length,@names) = mase($datafile); }	# MASE FORMAT 
else { ($seq_ref,$alignment_length,@names) = fasta($datafile); }				# FASTA FORMAT

%seq = %$seq_ref;

# TRANSLATE TO UPPERCASE TO AVOID PATTERN MATCHING ISSUES LATER AND GET RID OF UNUSUAL CHARACTERS
# 
print "\nChecking case. ";
my ($case, $converted);
foreach my $name (@names) {
	for (my $i=0; $i < @{$seq{$name}}; $i++) {
		if ($seq{$name}[$i] =~ /[acgt]/) { $seq{$name}[$i] =~ tr/acgt/ACGT/; $case++; }
		if ($seq{$name}[$i] !~ /[ACGT?-]/) { 
			unless (defined $converted) {
				$converted = "Converted $seq{$name}[$i] to ?.";
			}	
			$seq{$name}[$i] = "?"; 
		}
	}
}	
if (defined $case) { print "Converted lower-case bases to upper-case.\n"; }
else { print "Case is unchanged.\n"; }
if (defined $converted) { print "$converted\n"; }

## MERGE IDENTICAL SEQUENCES INTO A SINGLE HAPLOTYPE
#									# Note: Ref is undetermined when deciding haplotypes
if (exists $parameters{'M'}) {
	print "\n";
	my %seenhap;
	my @newnames;
	open HAPS, ">$hapfile" or die "couldn't open $hapfile : $!";
	open HAPSEQ, ">$hapseqfile" or die "couldn't open $hapseqfile : $!";
	foreach my $name1 (@names) {					# foreach seq
		if ($seenhap{$name1}) { next; }
		my @hapnames;
		foreach my $name2 (@names) {				# identify all matching seqs
			if ($name1 eq $name2) { next; }
			my $diff = "no";
			for (my $i=0; $i < @{$seq{$name1}}; $i++) {
			       if ($seq{$name1}[$i] eq $seq{$name2}[$i]) { next; }
			       elsif (($seq{$name1}[$i] eq "?") || ($seq{$name2}[$i] eq "?")) { next; }
			       elsif ($seq{$name1}[$i] ne $seq{$name2}[$i]) { $diff = "yes"; last; }
			}
			if ($diff eq "no") {
				push(@hapnames, $name2);
				$seenhap{$name2}++;	
			}
		}	       
		my $minamb = $alignment_length;				# identify which seq has the fewest ambiguous sites
		my $hapname;
		foreach my $name ($name1, @hapnames) {
			my $amb = 0;
			foreach (@{$seq{$name}}) { if ($_ eq "?") { $amb++; } }
		      	if ($amb < $minamb) { $hapname = $name; $minamb = $amb; }	
		}
		my @others;
		foreach (@hapnames) { if ($hapname ne "$_") { push (@others, $_); } }

		print "$hapname: @others\n";			
		print HAPS "$hapname: @others\n";

		push(@newnames, $hapname);
	}
	print "\nOriginal seqnames: @names\nNew merged hapnames: @newnames\n";
	print HAPS "\nOriginal seqnames: @names\nNew merged hapnames: @newnames\n";

	@names = @newnames;
}


## REFERENCE SEQUENCE(S) 
#						# $refseqnames are removed from @names index but %seq is kept
if (exists $parameters{'r'}) {			# > 1 refseq to be IDed using user provided $refpattern
	my @newnames;
	foreach (@names) { 
		if ($_ =~ /$refpattern/i) { push(@refseqnames, $_); }
		else { push (@newnames, $_); }
	}	
	@names = @newnames;
}
elsif (exists $parameters{'k'}) {		# create a majority consensus to use as reference
	push(@refseqnames, "consensus");
	for (my $i=0; $i < @{$seq{$names[0]}}; $i++) {
		my %seenbase;
		foreach my $name (@names) { $seenbase{$seq{$name}[$i]}++; }
		if (defined $seenbase{'?'}) {	# count ? as non-gaps
			foreach my $base (keys %seenbase) {
				if (($base eq "?") || ($base eq "-")) { next; }
				$seenbase{$base} += $seenbase{'?'};
			}
		}
		my $max = 0;
		my $majbase;
		foreach my $base (keys %seenbase) {		# check through all possibilities at each site in alignment in random order
			if ($seenbase{$base} > $max) { $majbase = $base; $max = $seenbase{$base}; }
		}
		push (@{$seq{"consensus"}},$majbase);
	}
}
else { push(@refseqnames, shift @names); }	# only 1st seq is treated as refseq

# print out merged haplotypes with majority consensus if these options have been selected
#
if (exists $parameters{'M'}) {
	my @printnames = @names;
	if (exists $parameters{'k'}) { @printnames = ("consensus", @names); }
	foreach (@printnames) {
		my $seq = join ("",@{$seq{$_}});
		print HAPSEQ ">$_\n$seq\n";
	}
	close HAPSEQ;
}

## CHECK REFERENCE SEQUENCE(S) FOR HOMOPOLYMER REPEATS 
#	
my ($previous, $base);
my $start = 0;
my $run = 1;
my $longrun = 0;
foreach my $refname (@refseqnames) {
	print "\nChecking for homopolymer runs (>= $hotspotlength) in $refname sequence (possible indel hotspots):\n";	

	for (my $i=0; $i < $alignment_length; $i++) {
		if ($seq{$refname}[$i] =~ /([acgt])/i) {
		       	unless (defined $previous) { $previous = $1; next; }	
			if ($1 eq $previous) { 
				if ($start == 0) { $start = $i; } 	# start of a run	
				$run++; 
				$base = $1;
			}
			else { 						# end of a run
				if ($run >= $hotspotlength) { 
					print "$datafile $refname homopolymer run ($base)n $run bp long at positions $start .. ".($i-1).".\n"; 
					$longrun++;
				}
				$previous = $1;				# rezero for next run search 
				$run = 1; 
				$start = 0;
			}

		}

	}
}
if ($longrun == 0) { print "$datafile @refseqnames has no homopolymer runs.\n"; }
	
## CREATE CODON POSITION INDEX FROM REFERENCE SEQUENCES
#
my @index;
my (%codon_pos,$started,$finished);
my ($conflicts, @conflicts, %seenconflicting, $conflictstart, $conflictend);

for (my $i=0; $i < $alignment_length; $i++) {
	my @with_seq;
	foreach my $refname (@refseqnames) {
		unless (defined $codon_pos{$refname}) { $codon_pos{$refname} = 0; }
		if ($seq{$refname}[$i] =~ /[?a-z]/i) { 			# ALL AMBIGUITY CODES ARE RECOGNISED FOR REFSEQ
			unless (defined $started) { $started = $i; }	# INDEX START
			$finished = $i;
			push(@with_seq,$refname); 
			$codon_pos{$refname} = codon_incr($codon_pos{$refname});
		}
	}
	unless (defined $started) { push(@index,"8"); }			# reference sequence has not started yet
	elsif (@with_seq == 0) { push(@index,"-"); }
	elsif (@with_seq == 1) { push(@index,$codon_pos{$with_seq[0]}); }
	elsif (@with_seq > 1) {						# > 1 reference sequence
		my $codon_pos;
		for (my $j=0; $j < @with_seq; $j++) {			# check if there's a conflict for the codon position
			if ($j==0) { $codon_pos = $codon_pos{$with_seq[0]}; }
			elsif ($codon_pos != $codon_pos{$with_seq[$j]}) {
				unless ($seenconflicting{$with_seq[0]}++) { push(@conflicts, $with_seq[0]); }
				unless ($seenconflicting{$with_seq[$j]}++) { push(@conflicts, $with_seq[$j]); }
				unless (defined $conflictstart) { $conflictstart = $i+1; }
				$conflictend = $i+1;
				$codon_pos = 0;				# call as 0 if there is disagreement on codon position and save details of refs involved
			}

		}
		push(@index,$codon_pos);
		if ($codon_pos == 0) { $conflicts++; }	
	}
}
if (defined $conflicts) { warn "\nReference sequences (@conflicts) show different codon positions at $conflicts positions ($conflictstart..$conflictend) : conflicts will be called as codon position 0.\n";	}

if (@index != $alignment_length) { die "problem making the index (".@index." bp) for the alignment ($alignment_length bp)\n"; }
if ($finished != $alignment_length - 1) { 				# reference sequence has finished
	for (my $i=$finished+1; $i < $alignment_length; $i++) {
		$index[$i] = 9;
	}
}

# PRINT A FASTA FILE WITH INDEX IF REQUESTED 
# 
if (exists $parameters{'f'}) {						

	open INDEXFA, ">$fastaindex" or die "couldn't open $fastaindex : $!";

	my $index = join ("",@index);
	my $nogaps;
	($nogaps = $index) =~ tr/-98//d;
	print INDEXFA ">index ".(length $nogaps)." bp\n$index\n";
	foreach my $name (@refseqnames,@names) {
		my $seq = join("",@{$seq{$name}});
		print INDEXFA ">$name\n$seq\n";
	}
}

print "\nNo. of reference sequences :\t".@refseqnames."\nReference sequences :\t\t@refseqnames\nNo. of sequences to analyse :\t".@names."\n";


## MAIN 
#
open OUT, ">$outfile" or die "couldn't open $outfile : $!";
open LENGTHS, ">$lengthfile" or die "couldn't open $lengthfile : $!";
if (exists $parameters{'I'}) { open INSERTIONS, ">$insertionfile" or die "couldn't open $insertionfile : $!"; } 


# HEADER
unless (exists $parameters{'h'}) { 
	print OUT "name\ttype\tvalue\tAstart .. Aend\tcodon_pos\n"; 
	print LENGTHS "name\tinitial_length\tfinal_length\tinsertion_length\tnon_unique_length\tfewseq_length\talphaI\talphaU\talphaU1st\talphaU2nd\talphaU3rd\talphaU0\talphaUa\talphaUc\talphaUg\talphaUt";

	if ($getuniqueatshared eq "yes") {
	#	print "\talphaU?";
		if (exists $parameters{'c'}) {
			foreach my $type (@neighbour_changes) { print LENGTHS "\talphaU$type"; }
		}
		else {
			foreach my $change (@changes,"?") { print LENGTHS "\talphaU$change"; }
		}
	}
	elsif (exists $parameters{'c'}) {
		foreach my $type (@neighbour_bases) { print LENGTHS "\talphaU$type"; }
	}
	print LENGTHS "\n";
}

my (%initial_length,%final_length,%insertion_length,%ins_allpos,%del_allpos,%ins_starts,%ins_ends,%ins_length,%start,%del_starts,%del_ends,%del_length,@allins_lengths,@alldel_lengths,@allcomplex_types,@allCdel_lengths,%ins_class,%del_class);
my (%unique_ps,%unique_ps_type,%unique_ps_codon,%shared_pos,%fewseq_length);
my (%ps_length,%alphaU,%alphaI,%alphaIJustin);
my (%complex_types,%complex_starts,%complex_ends,%Cdel_length,%complex_class);
my (%ambanc,%ambancestors);
my $unique_ps_count = 0;

foreach my $name (@names) {	

	$initial_length{$name} = 0;						# zero lengths at start to avoid nulls later
	$final_length{$name} = 0;
	$insertion_length{$name} = 0;
	$fewseq_length{$name} = 0;
	$ps_length{$name} = 0;
	$alphaU{$name} = 0;
	foreach my $codon_position (1,2,3,0) {
		$alphaU{"$name**$codon_position"} = 0;
	}
	foreach my $base ("A", "C", "G", "T", "?") {					# do not use qw with "?"
	       	$alphaU{"$name$base"} = 0;
	}
	foreach my $change (@changes) {
		$alphaU{"$name$change"} = 0;
	}
	if ((exists $parameters{'c'}) && ($getuniqueatshared eq "no")) {
		foreach my $type (@neighbour_bases) { $alphaU{"$name$type"} = 0; }
	}
	elsif (exists $parameters{'c'}) {
		foreach my $type (@neighbour_changes) { $alphaU{"$name$type"} = 0; }
	}

	for (my $i=0; $i < $alignment_length; $i++) {
			
		
# INSERTIONS DEFINED BY A "-" IN INDEX WHERE NOT AT UNDEFINED START OR END OF INDEX
# 
		if (($seq{$name}[$i] =~ /[a-z]/i) && ($index[$i] eq '-')) { 	# INSERTION
			push (@{$ins_allpos{$name}},$i);			# store every ins pos within orf region
		}								# index is 8 or 9 away from orf, never '-'

# DELETIONS DEFINED BY A "-" IN THIS SEQUENCE WHERE INDEX IS NOT "-" AND THIS SEQUENCE HAS STARTED
#
		if (!defined $start{$name}) {
			unless (exists $parameters{'t'}) { $start{$name} = 0; }	# COUNT 5' truncation AS DELS
				 						# DO NOT COUNT 5' truncation AS DELS
			elsif ($seq{$name}[$i] =~ /[a-z]/i) { $start{$name} = $i; }	
		}
										# DELETION
		if ((defined $start{$name}) && ($index[$i] ne '-') && ($seq{$name}[$i] eq '-')) {
			push (@{$del_allpos{$name}},$i);			# store every del pos
		}


# POINT SUBS
#
		if (($seq{$name}[$i] =~ /[a-z]/i) && ($index[$i] =~ /[0123]/)) {
			my (%bases,$not_unique,%point_sub, @others);
			my $count_other = 0;
			my $point_sub = 0;
			
			if ($includeref eq "yes") { @others = (@refseqnames, @names); }
			else { @others = @names; }
			
			foreach my $other (@others) {			
				if (($other ne $name) && ($seq{$other}[$i] =~ /[a-z]/i)) {
					$count_other++;
					$bases{$seq{$other}[$i]}++;		# grab other base states for ancestor

										# not unique flag
					if ($seq{$other}[$i] eq $seq{$name}[$i]) { $not_unique++; }
					else { 
						$point_sub++; 			# point sub flag
						$point_sub{$seq{$other}[$i]}++;	# are all other bases unique?
					}

				}
			}
			my %allbases = %bases;					# are all sites unique?
			$allbases{$seq{$name}[$i]}++;				# too little info to infer ancestry?
			my $not_all_unique;
			if (($count_other + 1) != (keys %allbases)) { $not_all_unique++; }
						
										# NEED AT LEAST 3 SEQS 
										# THAT ARE NOT ALL UNIQUE
			if (($count_other >= 2) && (defined $not_all_unique)) {	

				my $non_unique_flag;					
										# NON-UNIQUE POINT SUB
				if ((defined $not_unique) && ($point_sub > 1)) {
					my $not_all_others_unique;
					foreach my $base (keys %point_sub) {	
						if ($point_sub{$base} > 1) { $not_all_others_unique++; }
			#			print "hello: $base $not_all_others_unique\n";
					}
					if (defined $not_all_others_unique) {	# do not include site if all other bases 
									 	# are unique
						$shared_pos{$i}++;
						$non_unique_flag++;
					}
				}
										# UNIQUE POINT SUB
				elsif (($point_sub > 0) && (!defined $not_unique)) {	

					my @temp = ancestor(%bases);
					my $ancestral = shift @temp;

			#		print "hello1: $name [$i] $ancestral->$seq{$name}[$i]\n";
			
					push(@{$unique_ps{$name}}, $i);
					push(@{$unique_ps_codon{$name}}, $index[$i]); 
					if (exists $parameters{'c'}) {
						my ($nbefore,$nafter) = neighbours(\@{$seq{$name}},$i,\@index);
						
						push(@{$unique_ps_type{$name}}, "$nbefore$ancestral->$seq{$name}[$i]$nafter");
					}
					else { 
						push(@{$unique_ps_type{$name}}, "$ancestral->$seq{$name}[$i]"); 
					#	if ($ancestral =~ /\[/) { print "ERROR: $name\t$ancestral\t$i\n"; }

					}
					
					$unique_ps_count++;
				#	print OUT "$name\tunique_ps\t$ancestral->$seq{$name}[$i]\t$i .. $i\t$index[$i]\n";
				}

# POINT SUB LENGTH COUNTERS (SEQ FOR INDEX, FOR THIS STUDY SEQ AND AT LEAST 2 OTHER SEQS) SHARED DIFFS INCLUDED.
#	store lengths in %ps_length: key = nameindex_state length, or key = name for total_length
				
				$ps_length{"$name$index[$i]"}++;
				$ps_length{$name}++;

	# NUMBER OF AVAILABLE BASES FOR DIFFERENT TYPES OF UNIQUE POINT SUBSTITUTION
	# 			
										# GRAB ANCESTRAL STATE
				my ($ancestral, @non_uniquebases) = ancestor(%allbases);
				
										# EXCLUDING NON-UNIQUE POINT SUBS
				if ($getuniqueatshared eq "no") {
					unless (defined $non_unique_flag) {
										# do not count this base even if this
										# particular seq was unique.
						if (@non_uniquebases == 1) {


							$alphaU{"$name**$index[$i]"}++;	# codon position count
							$alphaU{$name}++;
							$alphaU{"$name$ancestral"}++;

										# AVAILABLE BASES with NEAREST NEIGHBOUR
							if (exists $parameters{'c'}) {
								my @neighbours = neighbours(\@{$seq{$name}},$i,\@index);
								my $at = atsummary(@neighbours);
								$alphaU{"$name$at$ancestral"}++;
							}
						}
						else { 
							$shared_pos{$i}++;	
						}

					}
				}
				else {						# NOT EXCLUDING NON-UNIQUE POINT SUBS
					my $seq = $seq{$name}[$i];
					
					$alphaU{$name}++;	
										# adjusted codon position count
					$alphaU{"$name**$index[$i]"} += (4 - @non_uniquebases)/3;

										# ONLY ONE NON-UNIQUE BASE. 
										# ANCESTOR UNAMBIGUOUS
					if (@non_uniquebases == 1) { 

						$alphaU{"$name$ancestral"}++; 

						foreach my $base (@bases) {
						
							if ($base =~ $ancestral) { next; }	# skip because not a change								
											# possible changes resolved
										# AVAILABLE BASES with NEAREST NEIGHBOUR
							if (exists $parameters{'c'}) {
								my @neighbours = neighbours(\@{$seq{$name}},$i,\@index);
								my $at = atsummary(@neighbours);
								$alphaU{"$name$at$ancestral->$base"}++;
							}
										# AVAILABLE BASES
							else { $alphaU{"$name$ancestral->$base"}++; }
						}							
					}

					else { 					# > 1 NON-UNIQUE BASE
						if (defined $non_unique_flag) {	# but base has not changed uniquely so do
										# not need to infer ancestor
							$alphaU{"$name$seq"}++;

							foreach my $base (@bases) {
											# could have been ancestral so 
								if ($base =~ /$ancestral/) { next; } # sub type not available								
										# possible changes resolved
										# AVAILABLE BASES with NEAREST NEIGHBOUR
								if (exists $parameters{'c'}) {
									my @neighbours = neighbours(\@{$seq{$name}},$i,\@index);
									my $at = atsummary(@neighbours);
									$alphaU{"$name$at$seq->$base"}++;
								}
										# AVAILABLE BASES
								else { $alphaU{"$name$seq->$base"}++; }
							}							
						}
						else {				# ANCESTOR AMBIGUOUS
							$alphaU{"$name?"}++; 
							if (exists $parameters{'c'}) {
								my @neighbours = neighbours(\@{$seq{$name}},$i,\@index);
								my $at = atsummary(@neighbours);
								$alphaU{"$name$at?"}++;
							}

										# gather data for sites that have
										# unique subs with ambiguous ancestors
										# for print (to assist manual phylogenetic)
										#  resolution
							push(@{$ambanc{$i}}, $name);
							$ambancestors{$i} = $ancestral;
						}
					}

				}
			}
			else { $fewseq_length{$name}++; }
		}


# LENGTHS. NEED COUNTERS TO DEAL WITH END EFFECTS OF DELETIONS THAT EXTEND BEYOND REFERENCE SEQUENCE.
# (length does not include extensions)
		if ((defined $start{$name}) 					# deals with truncation status (see DELETION)
			&& ($seq{$name}[$i] =~ /[a-z-]/i) && ($index[$i] =~ /[0123]/)) {
			$initial_length{$name}++;				# size at time of arrival
										# before deletions or insertions
		}
		if ((defined $start{$name})
			&& ($seq{$name}[$i] =~ /[a-z]/i) && ($index[$i] =~ /[0123-]/)) {
			$final_length{$name}++;					# size after deletions and insertions
		}		
		if ((defined $start{$name})
			&& ($seq{$name}[$i] =~ /[a-z]/i) && ($index[$i] =~ /-/)) {
			$insertion_length{$name}++;
		}
	}

# SUMMARISE SUBSTITUTIONS
#
										# SUMMARISE POINT SUBS
	if ($getuniqueatshared eq "no") {						# do not count unique ps at non-unique sites
		if (defined @{$unique_ps{$name}}) {
			my (@new_unique_ps,@new_unique_ps_type,@new_unique_ps_codon);
			for (my $i=0; $i < @{$unique_ps{$name}}; $i++) {
				if (exists $shared_pos{$unique_ps{$name}[$i]}) {
					#$alphaU{"$name$index[$unique_ps{$name}[$i]]"}--;
					#$alphaU{$name}--;
					$unique_ps_count--;
				}
				else {
					push (@new_unique_ps, $unique_ps{$name}[$i]);
					push (@new_unique_ps_type, $unique_ps_type{$name}[$i]);
					push (@new_unique_ps_codon, $unique_ps_codon{$name}[$i]);
				}
			}
			@{$unique_ps{$name}} = @new_unique_ps;
			@{$unique_ps_type{$name}} = @new_unique_ps_type;
			@{$unique_ps_codon{$name}} = @new_unique_ps_codon;
		}
	}

										# SUMMARISE INSERTIONS
	my ($ins_startsref,$ins_endsref,$ins_lengthref) = starts_and_length(\@{$ins_allpos{$name}},\@index,\@{$seq{$name}});
	@{$ins_starts{$name}} = @$ins_startsref;
	@{$ins_ends{$name}} = @$ins_endsref;
	@{$ins_length{$name}} = @$ins_lengthref;

										# SUMMARISE DELETIONS
	my ($del_startsref,$del_endsref,$del_lengthref) = starts_and_length(\@{$del_allpos{$name}},\@index,,\@{$seq{$name}});
	@{$del_starts{$name}} = @$del_startsref;
	@{$del_ends{$name}} = @$del_endsref;
	@{$del_length{$name}} = @$del_lengthref;

										# SUMMARISE COMPLEX EVENTS
	unless (exists $parameters{'C'}) {
		my (%remove_ins,@new_del_starts,@new_del_ends,@new_del_length,@new_ins_starts,@new_ins_ends,@new_ins_length);
		for (my $j = 0; $j < @{$del_starts{$name}}; $j++) {
			my $complex_event_flag;
			my $previous_pos = $del_starts{$name}[$j] - 1;
			my $next_pos = $del_ends{$name}[$j] + 1;
			for (my $k = 0; $k < @{$ins_ends{$name}}; $k++) {
				if ($ins_ends{$name}[$k] == $previous_pos) {	# insertion then deletion?
										# COMPLEX EVENT
										# store complex event characters
					push (@{$complex_types{$name}}, "+$ins_length{$name}[$k]-$del_length{$name}[$j]");
					push (@{$complex_starts{$name}}, $ins_starts{$name}[$k]);
					push (@{$complex_ends{$name}}, $del_ends{$name}[$j]);
					push (@{$Cdel_length{$name}}, $del_length{$name}[$j]);
					$remove_ins{$k}++;			# store which insertions need to be removed
					$complex_event_flag++;
										# PRINT INSERTION COMPONENT TO FA
					if ((exists $parameters{'I'}) && ($ins_length{$name}[$k] >= $insthreshold)) {
						my $ins_seq = join('',@{$seq{$name}}[$ins_starts{$name}[$k] .. $ins_ends{$name}[$k]]);
						print INSERTIONS ">$name $ins_length{$name}[$k] from complex event ";
						print INSERTIONS "-$del_length{$name}[$j]+$ins_length{$name}[$k] ";
						if ($ucscpos eq "no") { print INSERTIONS ($ins_starts{$name}[$k]+1)." .. ".($ins_ends{$name}[$k]+1)." in "; }
						else { print INSERTIONS "$ins_starts{$name}[$k] .. ".($ins_ends{$name}[$k]+1)." in "; }
						print INSERTIONS "$datafile\n$ins_seq\n";
					}
				}
				if ($ins_starts{$name}[$k] == $next_pos) {	# deletion then insertion?
										# COMPLEX EVENT
										# store complex event characters
					push (@{$complex_types{$name}}, "-$del_length{$name}[$j]+$ins_length{$name}[$k]");
					push (@{$complex_starts{$name}}, $del_starts{$name}[$j]);
					push (@{$complex_ends{$name}}, $ins_ends{$name}[$k]);
					push (@{$Cdel_length{$name}}, $del_length{$name}[$j]);
					$remove_ins{$k}++;			# store which insertions need to be removed
					$complex_event_flag++;
										# PRINT INSERTION COMPONENT TO FA
					if ((exists $parameters{'I'}) && ($ins_length{$name}[$k] >= $insthreshold)) {
						my $ins_seq = join('',@{$seq{$name}}[$ins_starts{$name}[$k] .. $ins_ends{$name}[$k]]);
						print INSERTIONS ">$name $ins_length{$name}[$k] from complex event ";
						print INSERTIONS "-$del_length{$name}[$j]+$ins_length{$name}[$k] ";
						if ($ucscpos eq "no") { print INSERTIONS ($ins_starts{$name}[$k]+1)." .. ".($ins_ends{$name}[$k]+1)." in "; }
						else  { print INSERTIONS "$ins_starts{$name}[$k] .. ".($ins_ends{$name}[$k]+1)." in "; }
						print INSERTIONS "$datafile\n$ins_seq\n";
					}

				}
			}
			unless (defined $complex_event_flag) { 
				push (@new_del_starts, $del_starts{$name}[$j]);
				push (@new_del_ends, $del_ends{$name}[$j]);
				push (@new_del_length, $del_length{$name}[$j]);
			}
		}
		for (my $i=0; $i < @{$ins_starts{$name}}; $i++) {
			unless (exists $remove_ins{$i}) { 
				push (@new_ins_starts, $ins_starts{$name}[$i]);
				push (@new_ins_ends, $ins_ends{$name}[$i]);
				push (@new_ins_length, $ins_length{$name}[$i]);
			}
		}
		
		@{$del_starts{$name}} = @new_del_starts;			# remove complex events from dels
		@{$del_ends{$name}} = @new_del_ends;
		@{$del_length{$name}} = @new_del_length;
		@{$ins_starts{$name}} = @new_ins_starts;			# remove complex events from ins
		@{$ins_ends{$name}} = @new_ins_ends;
		@{$ins_length{$name}} = @new_ins_length;
	}
	

	if ($alphaU{$name} == 0) { next; }					# skip sequences which do not align with ref (see also below)	

										# POOL ALL INSERTION LENGTHS FOR ALIGNMENT
	if (defined $ins_length{$name}[0]) { @allins_lengths = (@allins_lengths,@{$ins_length{$name}}); }

										# POOL ALL DELETION LENGTHS FOR ALIGNMENT
	if (defined $del_length{$name}[0]) { @alldel_lengths = (@alldel_lengths,@{$del_length{$name}}); }

										# POOL ALL COMPLEX TYPES & LENGTHS 
	if (defined $complex_types{$name}[0]) { 
		@allcomplex_types = (@allcomplex_types,@{$complex_types{$name}}); 
		@allCdel_lengths = (@allCdel_lengths,@{$Cdel_length{$name}});
	}



# PRINT LENGTHS TO LENGTH FILE
# 
	my @non_unique_seq = @{$seq{$name}}[(keys %shared_pos)];		# calculate length from summary
	my $non_unique_length = grep (/[a-z]/i, @non_unique_seq);
	
	print LENGTHS "$name\t$initial_length{$name}\t$final_length{$name}\t$insertion_length{$name}\t$non_unique_length\t$fewseq_length{$name}";
										# calculate length from summary
	$alphaI{$name} = ($final_length{$name} - $insertion_length{$name} + $initial_length{$name}) / 2;
	$alphaIJustin{$name} = ($final_length{$name} + $initial_length{$name}) / 2;			

	print LENGTHS "\t$alphaI{$name}\t$alphaU{$name}";
										# calculate length from summary
	my $test;
	if ($getuniqueatshared eq "no") {
		$test = $final_length{$name} - $insertion_length{$name} - $non_unique_length - $fewseq_length{$name};
	}
	else { $test = $final_length{$name} - $insertion_length{$name} - $fewseq_length{$name}; }
	if ($test ne $alphaU{$name}) {						# BUG CHECK
		print "WARNING. $name. length test alphaU on fly $alphaU{$name} not equal to length from summaries $test.\nCheck for bugs\n";
	}
	foreach my $codon_pos (1,2,3,0) {
		if (defined $ps_length{"$name$codon_pos"}) {
			print LENGTHS "\t".$alphaU{"$name**$codon_pos"};
		}
		else { print LENGTHS "\t0"; }

	}
	foreach my $base (@bases) {
		print LENGTHS "\t".$alphaU{"$name$base"};
	}

	if (($getuniqueatshared eq "no") && (exists $parameters{'c'})) {
		foreach my $type (@neighbour_bases) {
			print LENGTHS "\t".$alphaU{"$name$type"}; 
		}
	}
	
	if ($getuniqueatshared eq "yes") {
		if (exists $parameters{'c'}) {
			foreach my $type (@neighbour_changes) {
				print LENGTHS "\t".$alphaU{"$name$type"}; 
			}
		}
		else {
			foreach my $change (@changes,"?") {
				print LENGTHS "\t".$alphaU{"$name$change"};
			}
		}
	}	

	print LENGTHS "\n";

}

# SUMMARISE CODON POSITIONS OF NON_UNIQUE SITES
#
my $first = 0; my $second = 0; my $third = 0; my @other;
my @nonunique_pos = sort { $a <=> $b } (keys %shared_pos);
foreach my $codon_position (@index[@nonunique_pos]) {
	if ($codon_position == 1) { $first++; }
	elsif ($codon_position == 2) { $second++; }
	elsif ($codon_position == 3) { $third++; }
	else { push (@other, $codon_position); }
}

# IDENTIFY SHARED INSERTIONS, DELETIONS AND COMPLEX EVENTS

my $ins_class_ref = identify_shared(\%ins_starts,\%ins_ends,\%ins_length,"insertion");
my $del_class_ref = identify_shared(\%del_starts,\%del_ends,\%del_length,"deletion");
%ins_class = %$ins_class_ref;
%del_class = %$del_class_ref;
unless (exists $parameters{'C'}) {
	my $complex_class_ref = identify_shared(\%complex_starts,\%complex_ends,\%Cdel_length,"complex");
	%complex_class = %$complex_class_ref;
}

	
# PRINT SUBSTITUTIONS & count totals for unique codon positions
# 	
my $ufirst = 0; my $usecond = 0; my $uthird = 0; my @uother;
my $uins = 0; my $udels = 0; my $ucomplex = 0;

foreach my $name (@names) {

	if ($alphaU{$name} == 0) { 	# skip sequences which do not align with ref (see also above)	
		print "WARNING: Skipping $name because it has no sequence that aligns with the reference sequence.\n";
		next; 
	}


	for (my $i=0; $i < @{$ins_starts{$name}}; $i++) {			# PRINT INSERTIONS
		if ($ucscpos eq "no") { print OUT "$name\t$ins_class{$name}[$i]\t$ins_length{$name}[$i]\t".($ins_starts{$name}[$i]+1)." .. ".($ins_ends{$name}[$i]+1)."\tn/a\n"; }
		else { print OUT "$name\t$ins_class{$name}[$i]\t$ins_length{$name}[$i]\t$ins_starts{$name}[$i] .. ".($ins_ends{$name}[$i]+1)."\tn/a\n"; }

		if ($ins_class{$name}[$i] !~ /shared/) { $uins++; }
		
		if ((exists $parameters{'I'}) && ($ins_length{$name}[$i] >= $insthreshold)) {	# fasta file if requested
			my $ins_seq = join('',@{$seq{$name}}[$ins_starts{$name}[$i] .. $ins_ends{$name}[$i]]);
			if ($ucscpos eq "no") { print INSERTIONS ">$name $ins_length{$name}[$i] ".($ins_starts{$name}[$i]+1)." .. ".($ins_ends{$name}[$i]+1)." in $datafile\n$ins_seq\n"; }
			else { print INSERTIONS ">$name $ins_length{$name}[$i] $ins_starts{$name}[$i] .. ".($ins_ends{$name}[$i]+1)." in $datafile\n$ins_seq\n"; }

		}
	}
	
	for (my $i=0; $i < @{$del_starts{$name}}; $i++) {			# PRINT DELETIONS
		if ($ucscpos eq "no") { print OUT "$name\t$del_class{$name}[$i]\t$del_length{$name}[$i]\t".($del_starts{$name}[$i]+1)." .. ".($del_ends{$name}[$i]+1)."\tn/a\n"; }
		else  { print OUT "$name\t$del_class{$name}[$i]\t$del_length{$name}[$i]\t$del_starts{$name}[$i] .. ".($del_ends{$name}[$i]+1)."\tn/a\n"; }
		if ($del_class{$name}[$i] !~ /shared/) { $udels++; }
	}

	unless (exists $parameters{'C'}) {					# PRINT COMPLEX EVENTS
		if (defined @{$complex_starts{$name}}) {
			for (my $i=0; $i < @{$complex_starts{$name}}; $i++) {
				if ($ucscpos eq "no") { print OUT "$name\t$complex_class{$name}[$i]\t$complex_types{$name}[$i]\t".($complex_starts{$name}[$i]+1)." .. ".($complex_ends{$name}[$i]+1)."\tn/a\n"; }
				else  { print OUT "$name\t$complex_class{$name}[$i]\t$complex_types{$name}[$i]\t$complex_starts{$name}[$i] .. ".($complex_ends{$name}[$i]+1)."\tn/a\n"; }
				if ($complex_class{$name}[$i] !~ /shared/) { $ucomplex++; }

			}
		}
	}
											
	if (defined $unique_ps{$name}[0]) {					# PRINT UNIQUE POINT SUBS
		for (my $i=0; $i < @{$unique_ps{$name}}; $i++) {
			if ($ucscpos eq "no") { print OUT "$name\tunique_ps\t$unique_ps_type{$name}[$i]\t".($unique_ps{$name}[$i]+1)." .. ".($unique_ps{$name}[$i]+1)."\t$index[$unique_ps{$name}[$i]]\n"; }
			else  { print OUT "$name\tunique_ps\t$unique_ps_type{$name}[$i]\t$unique_ps{$name}[$i] .. ".($unique_ps{$name}[$i]+1)."\t$index[$unique_ps{$name}[$i]]\n"; }
			if ($index[$unique_ps{$name}[$i]] == 1) { $ufirst++; }
			elsif ($index[$unique_ps{$name}[$i]] == 2) { $usecond++; }
			elsif ($index[$unique_ps{$name}[$i]] == 3) { $uthird++; }
			else { push ( @uother, $index[$unique_ps{$name}[$i]]); }
		}
	}
}


# PRINT SUMMARY OF NON-UNIQUE POINT SUBS TO FILE
#
if ((exists $parameters{'s'}) && ((keys %shared_pos) > 0)) {	
	open NONUNIQUE, ">$nonuniquefile" or die "couldn't open $nonuniquefile : $!"; 
	print NONUNIQUE "# ".@nonunique_pos." non-unique positions, ".@names." sequences in $datafile alignment\n";
	print NONUNIQUE "# ";
	foreach (@nonunique_pos) { print NONUNIQUE ($_+1)." "; }	# using normal (base-1) coordinates for output fasta file, because output does not give a range so ucsc is less applicable
	print NONUNIQUE "\n"; 
	print NONUNIQUE "# @index[@nonunique_pos]\n";
	print NONUNIQUE "# codon position totals.\t1st: $first\t2nd: $second\t3rd: $third\tother: ";
	if (@other > 0) { print NONUNIQUE @other." (@other)\n"; }
	else { print NONUNIQUE "0\n"; }
	
	foreach my $name (@refseqnames,@names) {
		print NONUNIQUE ">$name\n";
		foreach my $pos (@nonunique_pos) {
			print NONUNIQUE "$seq{$name}[$pos]";	
		}
		print NONUNIQUE "\n";
	}
}

# PRINT SUMMARY OF SITES WHERE THERE WAS A UNIQUE SUB AND THE ANCESTOR WAS AMBIGUOUS
# 
unless (($getuniqueatshared eq "no") || (keys %ambanc == 0)) {
	open AMBANCESTOR, ">$ambancfile" or die "couldn't open $ambancfile : $!";

	my @ambanc_pos = sort { $a <=> $b } (keys %ambanc);

	print AMBANCESTOR "position\tNo.ofSeqsToCheck\tname\tbase\tpossible_ancestors\tSEE BELOW FOR A SUMMARY OF ALL DATA FOR JUST THESE SITES\n";
	foreach my $pos (@ambanc_pos){
		my @ambnames = @{$ambanc{$pos}};

		print AMBANCESTOR ($pos+1)."\t".@ambnames;	# using normal (base-1) coordinates for output fasta file, because output does not give a range so ucsc is less applicable
		foreach my $name (@ambnames) { print AMBANCESTOR "\t$name\t$seq{$name}[$pos]\t$ambancestors{$pos}"; }
		print AMBANCESTOR "\n";
	}
	
	print AMBANCESTOR "\nALL SEQS FOR SITES THAT NEED MANUAL CHECKING:\n";
	foreach my $pos (@ambanc_pos) {	print AMBANCESTOR "\t$pos"; }
	print AMBANCESTOR "\n";

	foreach my $name (@names) {
		print AMBANCESTOR "$name";
		foreach my $pos (@ambanc_pos) { print AMBANCESTOR "\t$seq{$name}[$pos]"; }
		print AMBANCESTOR "\n";	
	}
}


# SUMMARY TO STDOUT
#

if (@allins_lengths > 0) { 
	print "\nNo. of unique insertions :\t$uins\n"; 
	print "No. of all insertions :\t\t".@allins_lengths."\n";
	print "Lengths of all insertions :\t@allins_lengths\n"; 	
}
else { print "\nNo insertions.\n"; }
	
if (@alldel_lengths > 0) { 
	print "\nNo. of unique deletions :\t$udels\n"; 
	print "No. of all deletions :\t\t".@alldel_lengths."\n";
	my $complete = numeric(@alldel_lengths);
	print "No. of all complete deletions :\t$complete\n";
	print "Lengths of all deletions :\t@alldel_lengths\n"; 
}
else { print "\nNo deletions.\n"; }

unless (exists $parameters{'C'}) {
	if (@allcomplex_types > 0) { 
		print "\nNo. of unique complex events :\t\t$ucomplex\n"; 
		print "No. of all complex events :\t\t".@allcomplex_types."\n";
		my $complete = numeric(@allCdel_lengths);
		print "No. of all complete complex events :\t$complete\n";
		print "Lengths of all complex events :\t\t@allcomplex_types\n"; 

	}
	else { print "\nNo complex events.\n"; }
}


print "\nNon-unique point sub positions :\t".(keys %shared_pos)."\n";
print "Codon positions of non-unique sites :\t$first : $second : $third";
if (@other > 0) { print "\tothers: ".@other."\n"; }
else { print "\n"; }
print "No. of unique point subs :\t\t$unique_ps_count\n";
print "Codon positions of unique sites :\t$ufirst : $usecond : $uthird";
if (@other > 0) { print "\tothers: ".@uother."\n"; }
else { print "\n"; }
if (keys %ambanc > 0) { 
	print "Unique sub sites with >1 poss ancestor:\t".(keys %ambanc)." (check these & adjust all alphaUs).\n"; 
}
print "\n";



# CALCULATE RELATIVE RATES OF DIFFERENT TYPES OF EVENT
#
									# THIS HAS BEEN CHECKED BY RUNNING insertion THROUGH
									# relative_rate. GET THE SAME ANSWER AFTER REMOVING 
									# j SUBTRACTIONS & ADDITIONS
my ($lambdaI) = insertion(\@names,\%ins_length,\%unique_ps,\%alphaI,\%alphaU);
print "Relative rate of insertion (lambdaI) :\t$lambdaI ins per point sub\n";
my ($lambdaIJustin) = insertion(\@names,\%ins_length,\%unique_ps,\%alphaIJustin,\%alphaU);
print "Using Blumenstiel et al 02 :\t\t$lambdaIJustin ins per point sub\n\n";


my $maxdel_length = max(@alldel_lengths);				# GRAB MAX DEL LENGTH
									# ESTIMATE RELATIVE RATE OF DELS
my ($lambdaD,$lambdaDJustin) = relative_rate(\@names,$maxdel_length,\%del_length,\%unique_ps,\%final_length,\%initial_length,\%insertion_length,\%alphaU);
									# PRINT RELATIVE RATE RATE OF DELS
if ($lambdaD =~ /</) { print "No deletions. Probably < relative-rate given 1 event (1bp long).\n"; }
print "Relative rate of deletion (lambdaD) :\t$lambdaD dels per point sub\n";
print "Using Blumenstiel et al 02 :\t\t$lambdaDJustin dels per point sub\n\n";


unless (exists $parameters{'C'}) {

	my $maxCdel_length = max(@allCdel_lengths);			# GRAB COMPLEX EVENT MAX DEL LENGTH
									# ESTIMATE RELATIVE RATE OF COMPLEX EVENTS
	my ($lambdaC,$lambdaCJustin) = relative_rate(\@names,$maxCdel_length,\%Cdel_length,\%unique_ps,\%final_length,\%initial_length,\%insertion_length,\%alphaU);


	if ($lambdaC =~ /</) { print "No complex events. Probably < relative-rate given 1 event (1bp long).\n"; }
	print "Relative rate of complex events (lambdaC) :\t$lambdaC complex events per point sub\n";
	print "Using Blumenstiel et al 02 :\t\t\t$lambdaCJustin complex events per point sub\n";
}


					#################
					## SUBROUTINES ##
					#################	

				   #  IN ORDER OF APPEARANCE  #


## READ IN ALIGNMENT FROM NEXUS FORMAT 
#
sub nexus {

	my $datafile = shift;

	my (%seq, @names, %seen);
	my $length = 0;

	open DATA, "<$datafile" or die "couldn't open $datafile : $!";	

	while (<DATA>) {	
		if (/^(\S+)\s+([$seq_recog]+)\s*?\[*?\d*?\]*?$/mi) {
			my $name = $1;
			my $seq = $2;
			unless ($seen{$name}++) { push (@names, $name); }
			else { print "WARNING: $name has been seen $seen{$name} times. Only it's last sequence will be analysed.\n"; }
			if ($length == 0) { $length = length ($seq); }
			elsif (length($seq) != $length) { die "ERROR: not all sequences in $datafile are the same length: $1 is not $length bp\n"; }
			@{$seq{$name}} = split (//, $seq);
		}
	}
	print "Found ".@names." sequences in nexus format. All sequences are $length bp long.\n";

	return (\%seq,$length,@names);
}

## READ IN ALIGNMENT FROM MASE FORMAT 
#
sub mase {

	my $datafile = shift;
	my (%seq, @names, %seen, $name, $seq);
	my $length = 0;

	open DATA, "<$datafile" or die "couldn't open $datafile : $!";	

	my $getname = "no";
	my $getseq = "no";
	while (<DATA>) {	
		if (/^\;/m) {
			unless (/^\;\;/m) {
				if (defined $name) { 			# FINISH OFF PREVIOUS MASE ENTRY
					if ($length == 0) { $length = length($seq); }
					elsif (length($seq) != $length) { 
						die "ERROR: not all sequences in $datafile are the same length: $name is not $length bp\n$seq\n"; 
					}
					if ($seq =~ /[acgt]/i) { 
						unless ($seen{$name}++) { push (@names, $name); }
						else { print "WARNING: $name has been seen $seen{$name} times. Only its last sequence will be analysed.\n"; }
						@{$seq{$name}} = split (//, $seq); 	
						$seq = "";
					}
				}
				$getname = "yes"; 
				$getseq = "no";
			} 	 
			next;
		}
		elsif ($getname eq "yes") {
			/(\S+)/;
			$name = $1;
			$getname = "no";
			$getseq = "yes";
		}
		elsif ($getseq eq "yes") {
			/(\S+)/;
			$seq .= $1;	
		}

	}

	if ($seq =~ /[acgt]/i) { 				# FINISH OFF LAST MASE ENTRY
		unless ($seen{$name}++) { push (@names, $name); }
		else { print "WARNING: $name has been seen $seen{$name} times. Only its last sequence will be analysed.\n"; }
		@{$seq{$name}} = split (//, $seq); 	
	}
	print "Found ".@names." sequences in mase format. All sequences are $length bp long.\n";

	return (\%seq,$length,@names);

}

## READ IN ALIGNMENT FROM FASTA FORMAT 
#
sub fasta {

	my $datafile = shift;

	my (%seq, @names, $name, %seen, $seq);
	my $length = 0;

	open DATA, "<$datafile" or die "couldn't open $datafile : $!";	

	while (<DATA>) {
		if (/^>(\S+)/mi) {				# NAME LINE

			if (defined $name) { 			# FINISH OFF PREVIOUS FASTA ENTRY
				if ($length == 0) { $length = length($seq); }
				elsif (length($seq) != $length) { 
					die "ERROR: not all sequences in $datafile are the same length: $name is not $length bp\n"; 
				}
				if ($seq =~ /[acgt]/i) { 
					unless ($seen{$name}++) { push (@names, $name); }
					else { print "WARNING: $name has been seen $seen{$name} times. Only its last sequence will be analysed.\n"; }
					@{$seq{$name}} = split (//, $seq); 	
				}
				else { print "WARNING: No sequence for $name in fasta file. It is not included in the analysis.\n"; }
				$seq = '';
			}
			
			$name = $1;				# START NEW FASTA ENTRY
		}
		elsif (/^([$seq_recog]+)$/mi) { $seq .= $1; }	# GATHER SEQ (even if spread across multiple lines)
		elsif (/\S+/) { print "unrecognised line in datafile: [$_]\n"; }
	}
	if ($seq =~ /[acgt]/i) { 				# FINISH OFF LAST FASTA ENTRY
		if ($length == 0) { $length = length($seq); }
		elsif (length($seq) != $length) { 
			die "ERROR: not all sequences in $datafile are the same length: $name is not $length bp\n"; 
		}

		unless ($seen{$name}++) { push (@names, $name); }
		else { print "WARNING: $name has been seen $seen{$name} times. Only its last sequence will be analysed.\n"; }
		@{$seq{$name}} = split (//, $seq); 		
	}
	else { print "WARNING: No sequence for $name in fasta file. It is not included in the analysis.\n"; }


	print "Found ".@names." sequences in fasta format. All sequences are $length bp long.\n";

	return (\%seq,$length,@names);
}


# INCREMENT CODON POSITION (123123123etc)
# 
sub codon_incr {
	my $pos = shift;
	my $new_pos;

	if ((defined $pos) && ($pos == 3)) { $new_pos = 1; }
	else { $new_pos = ++$pos; }
	
	return $new_pos;
}


# SUMMARISE STARTS AND LENGTHS FROM AN ARRAY WITH EACH POSITION SHOWN 
# e.g. allpos: 2,3,4,5,36,40,41. starts: 2,36,40 lengths: 4,1,2 
sub starts_and_length {		

	my ($allpos_ref,$index_ref,$seq_ref) = @_;
	my @allpos = @$allpos_ref;
	my @index = @$index_ref;
	my @seq = @$seq_ref;
	my (@starts,@ends,@lengths,$previous_pos,$length);

	foreach my $pos (@allpos) {
		
		if (!defined $previous_pos) { 			# the first 
			push (@starts,$pos);			# new start and length
			$length = 1;
		}
		elsif ($pos == $previous_pos+1) { 		# continuing ..
			$length++; 
		}
		else {
								# skip if seq and index are both '-'
			my $intervening_index = join ('', @index[$previous_pos+1..$pos-1]);
			my $intervening_seq = join ('', @seq[$previous_pos+1..$pos-1]);
			unless ( ('-' x (length $intervening_index) eq $intervening_index) && ('-' x (length $intervening_seq) eq $intervening_seq) ) { 
				push (@lengths,$length);	# save final length for previous
				push (@ends,$previous_pos);	# save final end for previous
				push (@starts,$pos);		# new start and length
				$length = 1; 
			}
			else { $length++; }			# if the del is continuing .. so should it's length
		}
		$previous_pos = $pos; 

	}
	if (defined $previous_pos) {
		push (@lengths,$length);			# save final length for previous
		push (@ends,$previous_pos);			# save final end for previous
	}

# ALIGNMENT START AND END EFFECTS
# 
	for (my $i=0; $i < @starts; $i++) {
								# WHAT ABOUT INSERTIONS OUTSIDE ORF?
		unless (($index[$starts[$i]] =~ /[0123-]/) && ($index[$ends[$i]] =~ /[0123-]/)) {
			if (($starts[$i] == '0') ||			# start is before start of alignment
				($ends[$i] == (@index -1))) {		# end is after end of alignment			
				$lengths[$i] = ">$lengths[$i]";		# prepend with >
			}
			else { $lengths[$i] = "=$lengths[$i]"; }	# deletion is within the alignment but extends outside the ORF
		}
	
		if ($index[$starts[$i]] eq '9') {		# no indels can start after end of ref seq
			splice (@starts,$i,1);			# remove these from arrays
			splice (@ends,$i,1);
			splice (@lengths,$i,1);
			$i--;					# set counter back 1 because array has shrunk
		}
		elsif ($index[$ends[$i]] eq '8') {		# no indels can end before start of ref seq
			splice (@starts,$i,1);			# remove these from arrays
			splice (@ends,$i,1);
			splice (@lengths,$i,1);
			$i--;					# set counter back 1 because array has shrunk
		}
	}

	
# RETURN SUMMARY ARRAYS
# 
	return (\@starts,\@ends,\@lengths);
}

sub identify_shared {

	my %class;
	
	my ($starts_ref,$ends_ref,$lengths_ref,$default) = @_;
	my %starts = %$starts_ref;
	my %ends = %$ends_ref;
	my %lengths = %$lengths_ref;
	
	foreach my $name (keys %starts) {
		for (my $i=0; $i < @{$starts{$name}}; $i++) {
			
			my $shared;						# CHECK EVERY EVENT 

			unless (defined $starts{$name}[$i]) { next; }
			
			foreach my $othername (keys %starts) {

				if ($name eq $othername) { next; }

				
				for (my $j=0; $j < @{$starts{$name}}; $j++) {

					unless (defined $starts{$othername}[$j]) { next; }
			
					if ( ($starts{$name}[$i] == $starts{$othername}[$j])	
						&& ($ends{$name}[$i] == $ends{$othername}[$j]) ) {
										# numeric length	
						if (($lengths{$name}[$i] =~ /^\d+$/)
						&& ($lengths{$name}[$i] == $lengths{$othername}[$j])) {

							$shared++;		# FLAG EVENTS WHICH ARE SHARED
						}
										# non-numeric length
						elsif ($lengths{$name}[$i] eq $lengths{$othername}[$j]) {
							$shared++;		# FLAG EVENTS WHICH ARE SHARED
						}
					}
				}

			}
			if (defined $shared) { $class{$name}[$i] = "shared_$default"; }
			else { $class{$name}[$i] = $default; }
		}
	}
	
	return \%class;
}

sub neighbours {						

# FIND THE NEAREST NEIGHBOURS
# skips over neighbouring gaps if these are cause by insertions in other sequence. 
# Assumes the nearest base seen in the current sequence is the neighbour (not an inferred ancestor)
	
	my ($seqref,$i,$indexref) = @_;
	my @seq = @$seqref;
	my @index = @$indexref;
	my $nafter = "?";
	my $nbefore = "?";
	
	for (my $j=$i+1; $j < @seq; $j++) { 			# look in 3' direction
		if ($seq[$j] ne '-') { $nafter = $seq[$j]; last; }
	}

	for (my $j=$i-1; $j >= 0; $j--) { 			# look in 5' direction
		if ($seq[$j] ne '-') { $nbefore = $seq[$j]; last; }
	}
	
	return ($nbefore, $nafter);
}

sub atsummary {
	my @neighbours = @_;
	my $at = 0;
	
	foreach my $n (@neighbours) { 
		if ($n =~ /[at]/i) { $at++; } 
		elsif ($n !~ /[gc]/i) { $at = '?'; last; }
	}
	return $at;
}

sub ancestor {

	my %bases = @_;
	my (@non_uniquebases, $ancestral);			# only create ambiguous 
	foreach my $base (keys %bases) {			# ancestral for multiple non-unique
		if ($bases{$base} > 1) { push(@non_uniquebases, $base); }
	}
	if (@non_uniquebases == 1) { $ancestral = $non_uniquebases[0]; }
	else { my $temp = join(",",@non_uniquebases); $ancestral = "[$temp]"; }

	return ($ancestral, @non_uniquebases);
}



sub numeric {
	my @array = @_;
	my $n = grep (/^[0-9]+$/m, @array);
	return $n;
}


sub insertion {

	my ($namesref,$ins_lengthref,$unique_psref,$alphaIref,$alphaUref) = @_;
	my @names = @$namesref;	
	my %ins_length = %$ins_lengthref;
	my %unique_ps = %$unique_psref;
	my %alphaI = %$alphaIref;
	my %alphaU = %$alphaUref;
	my ($sum_nI, $sum_nU_x_alpha_ratio, $lambdaI, $sum_nU);		

	foreach my $name (@names) {
		my $nI = @{$ins_length{$name}};	
		my $nU = @{$unique_ps{$name}};			
		$sum_nI += $nI;		
		unless ($alphaU{$name} == 0) {
			$sum_nU_x_alpha_ratio += $nU * ($alphaI{$name}/$alphaU{$name});
		}
		
	#	print "$name $nI $nU $alphaI{$name} $alphaU{$name}\n";
		$sum_nU += $nU;
	}
	#print "\nsum_nI: $sum_nI sum_nU $sum_nU adjusted sum_nU $sum_nU_x_alpha_ratio\n";
								# note. there will be problems if sum nU is 0
	if ((!defined $sum_nU_x_alpha_ratio) || ($sum_nU_x_alpha_ratio == 0)) { 
		print "ERROR. No point subs in this dataset. Relative rates of insertion will fail.\n"; 
		$lambdaI = "error";
	}
	elsif ($sum_nI == 0) { 					# give a max for relative-rate if no insertions
		$sum_nI = 1;
		$lambdaI = "< ".($sum_nI / $sum_nU_x_alpha_ratio);
	}
	else { $lambdaI = $sum_nI / $sum_nU_x_alpha_ratio; }
	
	return ($lambdaI);
}

sub max {
	my @all = @_; 
	my $max;
	
	if (@all > 0) {						
		my @numeric = grep (/^[0-9]+/m, @all);
		@numeric = sort { $a <=> $b } @numeric;
		$max = $numeric[-1];
	}
	unless (defined $max) { $max = 0; }

	return $max;
}



sub relative_rate {

	my ($namesref,$maxj,$del_lengthref,$unique_psref,$final_lengthref,$initial_lengthref,$insertion_lengthref,$alphaUref) = @_;
	my @names = @$namesref;	
	my %del_length = %$del_lengthref;
	my %unique_ps = %$unique_psref;
	my %final_length = %$final_lengthref;
	my %initial_length = %$initial_lengthref;
	my %insertion_length = %$insertion_lengthref;
	my %alphaU = %$alphaUref;
	my ($no_events_flag, $lambda, $Justin_lambda, %sum_adjusted_nU, %Justin_sum_adjusted_nU);	

	if ($maxj == 0) { $maxj = 1; $no_events_flag++; }
	
	for (my $j = 1; $j <= $maxj; $j++) {
		my ($sum_n, $sum_nU);
		foreach my $name (@names) {
			my $n = grep (/^$j$/m, @{$del_length{$name}});
			my $nU = @{$unique_ps{$name}};
			my $alpha = (($final_length{$name} - $insertion_length{$name} + $j + $initial_length{$name})/2) - ($j + 1);
			my $Justin_alpha = (($final_length{$name} + $initial_length{$name})/2) - ($j + 1);
			$sum_n += $n;
			$sum_nU += $nU;

			unless ($alphaU{$name} == 0) {
				$sum_adjusted_nU{$j} += $nU * ($alpha / $alphaU{$name}); 
				$Justin_sum_adjusted_nU{$j} += $nU * ($Justin_alpha / $alphaU{$name}); 
			}
		}
		if (defined $no_events_flag) { $sum_n = 1; }
		if ((!defined $sum_adjusted_nU{$j}) || ($sum_adjusted_nU{$j} == 0)) { $lambda = "ERROR"; $Justin_lambda = "ERROR"; }
		else { 
			$lambda += $sum_n / $sum_adjusted_nU{$j};
			$Justin_lambda += $sum_n / $Justin_sum_adjusted_nU{$j};
		}
	#	print "del: $j sum_n: $sum_n sum_nU: $sum_nU sum_a_nU: $sum_adjusted_nU lambda_j: $lambda{$j}\n";
	}

	if (defined $no_events_flag) { $lambda = "< $lambda"; $Justin_lambda = "< $Justin_lambda"; }
	
	return ($lambda,$Justin_lambda);
}
