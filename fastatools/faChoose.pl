#!/usr/bin/perl


use warnings;
use strict;
use Getopt::Std;


my $program = 'faChoose.pl';				#name of script

my %parameters;						#input parameters
my ($infile,@names);
my $outfile = "faChoose.fasta";
my $seq_recog = "a-z0-9";				# recognises quality scores

getopts('i:n:o:e',\%parameters);

if (exists $parameters{"i"}) { $infile = $parameters{"i"}; }
if (exists $parameters{"n"}) { @names = split(/\s+/, $parameters{"n"}); }
if (exists $parameters{"o"}) { $outfile = $parameters{"o"}; }


unless ((exists $parameters{"i"}) && ($parameters{"n"})) {
	print "\n USAGE: $program -i <infile> -n <names>\n\n";

	print   "    -i\tfasta format infile\n";
	print   "    -n\tchosen name(s) of sequences e.g. 'S288c RM11'\n";
	print   "    -o\tname of outfile [$outfile]\n";
	print   "    -e\tallow exact matches only\n\n";
	exit;
}

my ($seq_ref,@namesinfile,%seq,$length);
($length,$seq_ref,@namesinfile) = fasta2strings($infile,$seq_recog);	# FASTA FORMAT SEQ
%seq = %$seq_ref;

open OUT, ">$outfile" or die "couldn't open $outfile : $!";
my @printed;
foreach my $fn (@namesinfile) {
	foreach my $un (@names) {	# print if name in the original file matches one of the names defined by the user
		if ((exists $parameters{'e'}) && ($fn eq "$un")) { print OUT ">$fn\n$seq{$fn}\n"; push (@printed,$fn); }
		elsif ((!defined $parameters{'e'}) && ($fn =~ /$un/mi)) { print OUT ">$fn\n$seq{$fn}\n"; push (@printed,$fn); }
	}
}

print "\nprinted ".@printed." sequences to $outfile\n@printed\n\n";
close OUT;

					#################
					## SUBROUTINES ##
					#################	

				   #  IN ORDER OF APPEARANCE  #


## READ IN ALIGNMENT FROM FASTA FORMAT 
#
# NOTE SEQ IS RETURNED AS A STRING NOT AS AN ARRAY, BETTER FOR LARGE DATA SETS
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
					$seq{$name} = $seq;	# STRINGS
					
					# @{$seq{$name}} = split (//, $seq); 	
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

		$seq{$name} = $seq;				# STRINGS

		# @{$seq{$name}} = split (//, $seq); 		
	}
	else { print "WARNING: No sequence for $name in fasta file. It is not included in the analysis.\n"; }


	print "\nFound ".@names." sequences in fasta format in $datafile. ";
	if (defined $alignment) { print "Sequences are different lengths\n\n"; }
	else { print "All sequences are $length bp long.\n\n"; }

	return ($length,\%seq,@names);
}


