#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

my $program = 'al2variable.pl';					#name of script

my $seq_recog = "0-9A-Z";					# SO INDEXES ARE RECOGNISED
my %parameters;							#input parameters
my $alignment;
my $prefix = "al2variable";

getopts('a:o:i:g',\%parameters);

if (exists $parameters{"a"}) { $alignment = $parameters{"a"}; }
if (exists $parameters{"o"}) { $prefix = $parameters{"o"}; }


unless (exists $parameters{"a"}) {
	print "\n USAGE: $program -a <gfa_file>\n\n";
	print   "    -a\tgapped fasta file of aligned seqs\n";
	print	"    -i\tstart pattern to match in names (eg 'T Q spar')[ALL]\n";
	print	"    -g\tinclude gaps\n";
	print 	"    -o\tfilename prefix for output [$prefix]\n\n";
	die;
}


## READ IN SEQUENCES FROM ALIGNMENT FILE 
#
my ($seq_ref,@names,%seq);
($seq_ref,@names) = fasta2strings($alignment,$seq_recog);	# FASTA FORMAT SEQ
%seq = %$seq_ref;


## GRAB INDEXES AND ALIGNMENT LENGTH
#
my ($rmindex,$index,@notindex,$alignmentlength);
foreach my $name (@names) { 
	if ($seq{$name} =~ /3/) { $index = $seq{$name}; print "\tindex = $name\n";  }
	elsif ($seq{$name} =~ /\d/) { $rmindex = $seq{$name}; print "\tRMindex = $name\n\n"; }
	else { 
		push (@notindex,$name); 
		$alignmentlength = length $seq{$name};
	}
}
@names = @notindex;


# IDENTIFY GROUP NAMES FROM THE COMMAND LINE ARGUMENTS
# 
my (@group, %group);
if (exists $parameters{'i'}) {
	my @patterns = split (/ /, $parameters{'i'});
	foreach my $pattern (@patterns) {
		foreach my $name (@names) {
			if ($name =~ /^$pattern/i) { push (@group, $name); $group{$name}++; }
		}
	}
	print "Study sequences: @group\n";
	@names = @group;
}

## PREPARE STUDY SEQUENCES
# 
my (%seqarray,%start,%end);
print "preparing sequence ... \n";
foreach my $name (@names) { 
	@{$seqarray{$name}} = split (//,$seq{$name}); 
	while ($seq{$name} =~ /[acgt]/igc) { 
		unless (defined $start{$name}) { $start{$name} = pos $seq{$name}; }
		$end{$name} = pos $seq{$name};	
	}
	print "$name\t$start{$name} .. $end{$name}\n";
}
print "done.\n\n";


## IDENTIFY VARIABLE SITES
# 
my (@variable, @ps, @gaps);
foreach (my $i=0; $i < $alignmentlength; $i++) {
	my $nt;
	foreach my $name (@names) {
		if (($i < ($start{$name}-1)) || ($i >= $end{$name})) { next; }	# skip start and end with no seq
		if ($seqarray{$name}[$i] =~ /[acgt]/i) {			# grab variable sites not including indels
			unless (defined $nt) { $nt = $seqarray{$name}[$i]; }
			elsif ($seqarray{$name}[$i] !~ /$nt/i) { 
				push (@variable,$i); 
				if ($nt eq '-') { push (@gaps,$i); }
				else { push (@ps,$i); } 
				last; 
			}
		}
										# grab variable sites with indels
		elsif ((exists $parameters{'g'}) && ($seqarray{$name}[$i] =~ /-/i)) {	 	
			unless (defined $nt) { $nt = $seqarray{$name}[$i]; }
			elsif ($seqarray{$name}[$i] !~ /$nt/i) { push (@variable,$i); push (@gaps,$i); last; }
		}

	}
}
print "variable sites:\t".@variable."\n";
print "point subs:\t".@ps."\n";
if (exists $parameters{'g'}) { print "gaps:\t\t".@gaps."\n"; }
else { print "Gaps ignored\n"; }


## CREATE VARIABLE SITE SEQS
#
my %varseq;
open OUT, ">$prefix.list" or die "couldn't open $prefix.list :$!";
print OUT "newpos\toldpos\n";
my $count = 0;
foreach my $i (@variable) {
	$count++;
	print OUT "$count\t".($i+1)."\n";
	foreach my $name (@names) { $varseq{$name} .= $seqarray{$name}[$i]; }
}
close OUT;

## PRINT FASTA FILE OF ONLY VARIABLE SITES
#
open OUT, ">$prefix.gfa" or die "couldn't open $prefix.gfa : $!";
foreach my $name (@names) {
	print OUT ">$name\n$varseq{$name}\n";
}
close OUT;


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

