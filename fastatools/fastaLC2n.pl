#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

my $program = 'fastaLC2n.pl';					#name of script

my $seq_recog = "0-9A-Z ";					# SO QUALITY FILES ARE RECOGNISED
my %parameters;							#input parameters
my ($file, $prefix,$outfile);

getopts('i:o:NA',\%parameters);


if (exists($parameters{"i"})) { $file = $parameters{"i"}; }
if (exists($parameters{"o"})) { $outfile = $parameters{"o"}; }


unless ((exists $parameters{"i"}) && (exists $parameters{"o"})) {
	print "\n $program converts lowercase sequence to 'n's\n";
	print "\n USAGE: $program -i <fasta_file> -o <outfile>\n\n";
	print   "    -i\tfasta file of one or multiple sequences or alignment scores\n";
	print 	"    -N\tconvert to uppercase N\n";
	print 	"    -A\tconvert ambiguity codes [KMRSWYVHDBkmrswyvhdb] to uppercase N\n";
	print 	"    -o\tname for output file\n\n";
	exit;
}


## READ IN SEQUENCES FROM ALIGNMENT FILE 
#
my ($seq_ref,@names,%seq);
($seq_ref,@names) = fasta2strings($file,$seq_recog);	# FASTA FORMAT SEQ
%seq = %$seq_ref;

open OUT, ">$outfile" or die "couldn't open $outfile : $!"; 
foreach my $name (@names) {
	if (exists $parameters{'A'}) { 	$seq{$name} =~ s/[KMRSWYVHDB]/N/g; }
	if (exists $parameters{'N'}) { 	$seq{$name} =~ s/[a-z]/N/g; }
	else { $seq{$name} =~ s/[a-z]/n/g; }
	print OUT ">$name\n$seq{$name}\n";

}
if (exists $parameters{'A'}) { 	print "Converted any [KMRSWYVHDB] to N\n"; }
if (exists $parameters{'N'}) { 	print "Converted any [a-z] to N\n"; }


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

	return (\%seq,@names);
}

