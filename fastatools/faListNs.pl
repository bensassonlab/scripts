#!/usr/bin/perl


use warnings;
use strict;
use Getopt::Std;


my $program = 'faFindNs.pl';				#name of script

my %parameters;						#input parameters
my ($infile);
my $seq_recog = "a-z0-9";				# recognises quality scores

getopts('i:',\%parameters);

if (exists $parameters{"i"}) { $infile = $parameters{"i"}; }

unless (exists $parameters{"i"}) {
	print "\n USAGE: $program -i <infile>\n\n";

	print   "    -i\tfasta format infile\n\n";
	exit;
}

my ($seq_ref,@namesinfile,%seq,$length);
($length,$seq_ref,@namesinfile) = fasta2strings($infile,$seq_recog);	# FASTA FORMAT SEQ
%seq = %$seq_ref;

foreach my $name (@namesinfile) {
#	print "Searching for \"N\" on $name ..\n";
	while ($seq{$name} =~ /(N+)/g) {
		my $nlength = length $1;
		my $end = pos ($seq{$name}) + 1;
		my $start = $end - $nlength;
		print "$name N_length: $nlength N_pos: $start..$end\n"; 
	}
}
print "\n";

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


