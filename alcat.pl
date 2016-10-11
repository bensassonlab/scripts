#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

my $program = 'alcat.pl';					#name of script

my $seq_recog = "0-9A-Z";					# SO INDEXES ARE RECOGNISED
my %parameters;							#input parameters
my ($infiles);
my $outfile = "alcat.gfa";

getopts('i:o:r:',\%parameters);

if (exists $parameters{"i"}) { $infiles = $parameters{"i"}; }
if (exists $parameters{"o"}) { $outfile = $parameters{"o"}; }
#if (exists $parameters{"r"}) { $refname = $parameters{"r"}; }

unless (exists $parameters{"i"}) {
	print "\n USAGE: $program -i '<fastafile> <fastafile> .. etc'\n\n";
	print   "    -i\tgfa alignments with the same taxa\n";
#	print	"    -r\tname of ref sequence for maximum length\n";
	print 	"    -o\tfilename for output [$outfile]\n\n";
	die;
}

my @infiles = split(/\s+/,$infiles);

# READ IN ALL SEQUENCE NAMES FROM ALL FILES 
#
my %allnames;
foreach my $infile (@infiles) {
	open DATA, "<$infile" or die "couldn't open $infile : $!";	
	while (<DATA>) { if (/^>(\S+)/mi) { $allnames{$1}++; } }
}
					
## READ IN SEQUENCES FROM ALIGNMENT FILE 
#
my (%catseq,$seq_ref,@names,%seq,$length,%seen,%names);

foreach my $infile (@infiles) {
	my @missing;
	($seq_ref,$length,@names) = fasta2strings($infile,$seq_recog);	# FASTA FORMAT SEQ
	%seq = %$seq_ref;
	@{$names{$infile}} = @names;

	foreach my $name (sort keys %allnames) {
		$seen{$name}++;
		unless (defined $seq{$name}) { 
			$catseq{$name} .= "N" x $length;
			push (@missing, $name);
		}
		else { 
			if (length($seq{$name}) < $length) { $seq{$name} .= "N" x ($length - length($seq{$name})); }
			$catseq{$name} .= $seq{$name}; 
		}
	}
	if (defined $missing[0]) { print "Missing from $infile:\t@missing\n"; }
}

## GRAB AND PRINT NEW FILE
#
open OUT, ">$outfile" or die "couldn't open $outfile : $!\n";
foreach my $name (sort keys %allnames) {
	print OUT ">$name\n$catseq{$name}\n";
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
					if (length($seq) > $length) { $length = length($seq); } # report the max seq length 

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

		if (length($seq) != $length) { 
			$alignment = 'no';
			if (length($seq) > $length) { $length = length($seq); } # report the max seq length 

		}

		$seq{$name} = $seq;					# STRINGS

		#@{$seq{$name}} = split (//, $seq); 			# ARRAYS
	}
	else { print "WARNING: No sequence for $name in fasta file. It is not included in the analysis.\n"; }


	print "\n$datafile :\t".@names." seqs. ";
	if (defined $alignment) { 
		warn "Sequences were different lengths in $datafile, so filled in the ends with Ns\n\n"; 
	}
	else { print "length: $length bp.\n"; }

	return (\%seq,$length,@names);
}



