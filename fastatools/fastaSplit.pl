#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

my $program = 'fastaSplit.pl';					#name of script

my $seq_recog = "0-9A-Z ";					# SO QUALITY FILES ARE RECOGNISED
my %parameters;							#input parameters
my $file;

getopts('i:',\%parameters);


if (exists($parameters{"i"})) { $file = $parameters{"i"}; }


unless (exists $parameters{"i"}) {
	print "\n $program splits a multi-fasta file into multiple fasta files\n\n";
	print "\n USAGE: $program -i <fasta_file>\n\n";
	print   "    -i\tfasta file of multiple sequences\n";
	print 	"\n";
	exit;
}


## READ IN SEQUENCES FROM ALIGNMENT FILE 
#
my ($seq_ref,@names,%seq);
($seq_ref,@names) = fasta2strings($file,$seq_recog);	# FASTA FORMAT SEQ
%seq = %$seq_ref;

foreach my $name (@names) {
	my $outname = $name;
	if ($name =~ /(\S+?)[:]/) { $outname = $1; }	# truncate name at special characters
	open OUT, ">$outname.fa" or die "couldn't open $name.fa : $!"; 
	print OUT ">$outname\n$seq{$name}\n";
	close OUT;

}


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

