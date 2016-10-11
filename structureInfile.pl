#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

my $program = 'structureInfile.pl';				#name of script

my $seq_recog = "0-9A-Z";					# SO INDEXES ARE RECOGNISED
my %parameters;							#input parameters
my ($alignment,$suffix,$popfile);
my $prefix = "al2variable";

my %bases = (
	a => 1,
	c => 2,
	g => 3,
	t => 4,
	n => -9,
	'-' => -9,
);


getopts('a:s:o:i:p:g',\%parameters);

if (exists $parameters{"a"}) { $alignment = $parameters{"a"}; }
if (exists $parameters{"s"}) { $suffix = $parameters{"s"}; }
if (exists $parameters{"o"}) { $prefix = $parameters{"o"}; }
if (exists $parameters{"p"}) { $popfile = $parameters{"p"}; }


unless ((exists $parameters{"a"}) || ($parameters{"s"})) {
	print "\n USAGE: $program -a <gfa_file>\n";
	print "    OR: $program -s <suffix>\n\n";

	print   "    -a\tgapped fasta file of aligned seqs\n";
	print   "    -s\tsuffix for multiple gapped fasta files e.g. gfa\n";
	print	"    -g\tinclude gaps\n";
	print	"    -p\tfile with list of popns: 'popNo StartOfName'\n";
	print 	"    -o\tfilename prefix for output [$prefix]\n\n";
	exit;
}

my (@alignments, @allnames, %varseq);


## READ IN ALIGNMENT FILES FROM CURRENT DIRECTORY 
#
if (defined $suffix) { 				# read directory of multiple gfa files

	opendir(DIR, ".") or die "can't opendir . : $!";

	while (defined(my $file = readdir(DIR))) {				# read each alignment file
		unless ($file =~ /(\S+?)\.$suffix$/m) { next; }
		my $name = $1;

#		print "Found $file : $name .. \n";
		push (@alignments, $file);
	}

}
else { push (@alignments, $parameters{'a'}); }


#### ALL FILES ARE USED IN THE STRUCTURE INFILE
#
#
#
#
#
my (%variable,%longseq,%alignmentlength);

foreach my $alignment (@alignments) {

	## READ IN SEQUENCES FROM ALIGNMENT FILE 
	#
	my ($seq_ref,@names,%seq);
	($seq_ref,@names) = fasta2strings($alignment,$seq_recog);	# FASTA FORMAT SEQ
	%seq = %$seq_ref;
	unless (@allnames) { @allnames = @names; }
	elsif (@names ne @allnames) { 
		die "ERROR: the ".@names." names in $alignment (@names) do not match the ".@allnames." names seen in other files (@allnames).\n"; 
	}

	## ALIGNMENT LENGTH AND CONVERT TO LOWER CASE
	#
	my $alignmentlength;
	foreach my $name (@names) { 
		$alignmentlength = length $seq{$name};
		$seq{$name} =~ tr/[A-Z]/[a-z]/;			# convert all to lowercase
		$alignmentlength{$alignment} = $alignmentlength;
	}

	## PREPARE STUDY SEQUENCES
	# 
	my (%seqarray,%start,%end);
#	print "preparing sequence ... \n";
	foreach my $name (@names) { 
		@{$seqarray{$name}}  = split (//,$seq{$name}); 	# add the seq from this alignment to any pre-existing seq
		push (@{$longseq{$name}}, @{$seqarray{$name}});
		while ($seq{$name} =~ /[acgt]/igc) { 
			unless (defined $start{$name}) { $start{$name} = pos $seq{$name}; }
			$end{$name} = pos $seq{$name};	
		}
		unless (defined $start{$name}) { $start{$name} = 0; $end{$name} = $alignmentlength; }	# where there is only missing sequence
#		print "$name\t$start{$name} .. $end{$name}\n";
	}
#	print "done.\n\n";


	## IDENTIFY VARIABLE SITES
	# 
	my (@ps, @gaps);
	foreach (my $i=0; $i < $alignmentlength; $i++) {
		my ($nt, $first);
		foreach my $name (@names) {
			if (($i < ($start{$name}-1)) || ($i >= $end{$name})) { next; }	# skip start and end with no seq
			if ($seqarray{$name}[$i] =~ /[acgt]/i) {			# grab variable sites not including indels
				unless (defined $nt) { $nt = $seqarray{$name}[$i]; }	# first sequence
				elsif ($seqarray{$name}[$i] !~ /$nt/i) { 		# nucleotide is not the same as in prev OTU
					push (@{$variable{$alignment}},$i); 
					if ($nt eq '-') { push (@gaps,$i); }
					else { push (@ps,$i); } 
					last; 
				}
			}
										# grab variable sites with indels
			elsif ((exists $parameters{'g'}) && ($seqarray{$name}[$i] =~ /-/i)) {	 	
				unless (defined $nt) { $nt = $seqarray{$name}[$i]; }
				elsif ($seqarray{$name}[$i] !~ /$nt/i) { push (@{$variable{$alignment}},$i); push (@gaps,$i); last; }
			}

		}
	}
	if (!defined $variable{$alignment}[0]) {
		print "There are no variable sites in $alignment. Skipping ..\n"; next;
	}
	print "variable sites:\t".@{$variable{$alignment}}."\n";
	print "point subs:\t".@ps."\n";
	if (exists $parameters{'g'}) { print "gaps:\t\t".@gaps."\n"; }
	else { print "Gaps ignored\n"; }


	## CREATE VARIABLE SITE SEQS
	#
	open OUT, ">$prefix.list" or die "couldn't open $prefix.list :$!";
	print OUT "newpos\toldpos\n";
	my $count = 0;
	foreach my $i (@{$variable{$alignment}}) {
		$count++;
		print OUT "$count\t".($i+1)."\n";
		foreach my $name (@names) { $varseq{$name} .= $seqarray{$name}[$i]; }
	}
	close OUT;

}

## READ FILE DESCRIBING POPNS IF ONE IS SPECIFIED
#
my %pops;
my $maxpop = 0;
if (defined $popfile) {
	open POP, "<$popfile" or die "couldn't open $popfile : $!";
	while (<POP>) {
		if (/(\d+)\s+(\S+)/) { 
			if ($1 > $maxpop) { $maxpop = $1; }		
			$pops{$2} = $1;	
		}
		else { print "warning: unrecognised format in $popfile. Will not use population information.\n"; }
	}
	$maxpop++;
}	


## PRINT FASTA FILE OF ONLY VARIABLE SITES
#
open OUT, ">$prefix.fasta" or die "couldn't open $prefix.fasta : $!";
foreach my $name (@allnames) {
	print OUT ">$name\n$varseq{$name}\n";
}
close OUT;



## PRINT STRUCTURE INFILE OF ONLY VARIABLE SITES
#
open INFILE, ">$prefix.infile" or die "couldn't open $prefix.infile : $!";

print INFILE "\t\t-1";	# print distance between loci
print "\n";
my $previous;

my $vscount = 0;
foreach my $alignment (@alignments) {

	my $first;
	$vscount += @{$variable{$alignment}};
	print "$alignment : ".@{$variable{$alignment}}." variable sites\n"; 
	foreach my $i (@{$variable{$alignment}}) {
		my $d;
		if (!defined $previous) { $previous = $i+1; next; }
		unless ($first++) { $d = -1; } 	
		else { $d = $i+1 - $previous; }
	#	print "$d = $i +1\t$previous\n";
		print INFILE "\t$d";
		$previous = $i;
	}

}
print INFILE "\n";

if (defined $popfile) { print "\nPopulations:\n\tName\tPop\tFlag\n"; }

foreach my $name (@allnames) {
	print INFILE "$name";			# print strain and genotypes
	if (defined $popfile) {			# Print Pop and Flag if popfile is provided
		my $popinfo = "no";
		foreach my $prename (keys %pops) {
			if ($name =~ /$prename/) { 
				print INFILE "\t$pops{$prename}\t1"; 
			       	print "\t$name\t$pops{$prename}\t1\n";
				$popinfo = "yes";
				last; 
			}	
		}
		if ($popinfo eq 'no') {
		       	print INFILE "\t$maxpop\t0";	
			print "\t$name\t$maxpop\t0\n";
		}		
	}
	else { print INFILE "\t8"; }		# Otherwise, add a dummy population
	my $length = 0;

	foreach my $alignment (@alignments) {
		foreach my $i (@{$variable{$alignment}}) {
			my $base = $longseq{$name}[$i+$length];
			my @test = keys %bases;
			unless (exists $bases{$base}) { 
				die "ERROR: unrecognised base for $name in $alignment '$base', @test. This is at position ".($i+$length)." in ".@{$longseq{$name}}." bp long sequence (length: $length)\n"; 
			}
			print INFILE "\t$bases{$base}";
		}
		$length += $alignmentlength{$alignment};
	}
	print INFILE "\n";

}

close INFILE;

print "\nstructure infile with ".@allnames." taxa and $vscount variable sites printed to $prefix.infile\n\n";

					#################
					## SUBROUTINES ##
					#################	

				   #  IN ORDER OF APPEARANCE  #


## READ IN ALIGNMENT FROM FASTA FORMAT 
#
# NOTE SEQ IS RETURNED AS AN ARRAY NOT AS A STRING, SLOW - FOR DETAILED ANALYSIS
sub fasta2strings {

	my $datafile = shift;
	my $seq_recog = shift;

	my (%seq, @names, $name, %seen, $seq, $alignment);
	my $length = 0;

	open DATA, "<$datafile" or die "couldn't open $datafile : $!";	

	while (<DATA>) {
		if (/^>(\S+)/mi) {				# NAME LINE

			if (defined $name) { 			# FINISH OFF PREVIOUS FASTA ENTRY
				if ($length == 0) { $length = length($seq); }
				elsif (length($seq) != $length) { 
					$alignment = 'no';
				}
				if ($seq =~ /[$seq_recog]/i) { 
					unless ($seen{$name}++) { push (@names, $name); }
					else { print "WARNING: $name has been seen $seen{$name} times. Only its last sequence will be analysed.\n"; }
					$seq{$name} = $seq;			# STRINGS
					
					#@{$seq{$name}} = split (//, $seq); 	# ARRAYS
				}
				else { print "WARNING: No sequence for $name in fasta file. It is not included in the analysis.\n"; }
				$seq = '';
			}
			
			$name = $1;				# START NEW FASTA ENTRY
		}
		elsif (/^([$seq_recog\-]+)\s*$/mi) { $seq .= $1; }	# GATHER SEQ (even if spread across multiple lines)
		elsif (/\S+/) { print "unrecognised line in datafile: [$_]\n"; }
	}
	if ($seq =~ /[$seq_recog]/i) { 				# FINISH OFF LAST FASTA ENTRY
		unless ($seen{$name}++) { push (@names, $name); }
		else { print "WARNING: $name has been seen $seen{$name} times. Only its last sequence will be analysed.\n"; }

		$seq{$name} = $seq;					# STRINGS

		#@{$seq{$name}} = split (//, $seq); 			# ARRAYS
	}
	else { print "WARNING: No sequence for $name in fasta file. It is not included in the analysis.\n"; }


	print "\n$datafile :\t".@names." seqs. ";
	if (defined $alignment) { die "Sequences are different lengths in $datafile\n\n"; }
	else { print "length: $length bp.\n"; }

	return (\%seq,@names);
}

