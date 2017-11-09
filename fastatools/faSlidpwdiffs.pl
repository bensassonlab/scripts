#!/usr/bin/perl


use warnings;
use strict;
use Getopt::Std;


my $program = 'faSlidpwsdiffs.pl';				#name of script

my %parameters;						#input parameters
my ($infile,$ref);
my $outfile;
my $seq_recog = "a-z";				
my $W = 5000;						# window size

my %ambcodes = (
	"M" => "AC",
	"R" => "AG",
	"W" => "AT",
	"S" => "CG",
	"Y" => "CT",
	"K" => "GT",
);

getopts('i:r:W:a',\%parameters);

if (exists $parameters{"i"}) { $infile = $parameters{"i"}; }
if (exists $parameters{"r"}) { $ref = $parameters{"r"}; }
if (exists $parameters{"W"}) { $W = $parameters{"W"}; }


unless ((exists $parameters{"i"}) && ($parameters{"r"})) {
	print "\n USAGE: $program -i <infile> -r <ref>\n\n";

	print   "    -i\tfasta format infile\n";
	print   "    -r\treference seq name (e.g. P34048)\n";
#	print   "    -a\tambiguity codes treated as both bases (Y=C,Y=T)\n";
	print   "    -W\tsliding window size [$W bp]\n\n";
	exit;
}

$outfile = "$ref$infile.tsv";

my ($seq_ref,@namesinfile,%seq,$length);
($length,$seq_ref,@namesinfile) = fasta2strings($infile,$seq_recog);	# FASTA FORMAT SEQ
%seq = %$seq_ref;

unless (defined $seq{$ref}) { die "error: could not find reference sequence called $ref\n"; }

open OUT, ">$outfile" or die "couldn't open $outfile : $!";
print OUT "file\tref\tstrain\tpos\tdiffs\tlengthNotN\tpDiff\n";
								 
foreach my $name (@namesinfile) { 
	if ($name eq $ref) { next; }					# skip the ref (compare other seqs to ref only)
	my $Wi = 0; my $diff = 0; my $T = 0; my $l = 0;				# position in window
	for (my $i=0; $i < @{$seq{$name}}; $i++) {
		if (($seq{$name}[$i] eq "N") || ($seq{$name}[$i] eq "n") || ($seq{$ref}[$i] eq "N") || ($seq{$ref}[$i] eq "n")) { $Wi++; }	# skip low quality sequence in length
		elsif ($seq{$name}[$i] eq $seq{$ref}[$i]) { $Wi++; $l++; }	# no difference from reference
		elsif (($seq{$name}[$i] =~ /[MRWSYK]/) && ($seq{$ref}[$i] =~ /[$ambcodes{$seq{$name}[$i]}]/)) {	# seq is heterozygous and matches one allele of ref 
		#	print "seq: $seq{$name}[$i]\t$seq{$ref}[$i] = not a difference\n";
		}
		elsif (($seq{$ref}[$i] =~ /[MRWSYK]/) && ($seq{$name}[$i] =~ /[$ambcodes{$seq{$ref}[$i]}]/)) {	# seq is heterozygous and matches one allele of ref 
		#	print "seq: $seq{$name}[$i]\t$seq{$ref}[$i] = not a difference\n";
		}
		else { $diff++; $Wi++; $l++; }
		if ($Wi	== $W) {					# reached the end of a window	
			$T += $Wi;					# keep track of the true position
			if ($l == 0) { print OUT  "$infile\t$ref\t$name\t$T\t$diff\t$l\tNA\n"; }
			else {	
				my $pDiff = $diff/$l;
				print OUT "$infile\t$ref\t$name\t$T\t$diff\t$l\t$pDiff\n";
			}
			$diff = 0; $Wi = 0; $l = 0;
		}
	}	
}
close OUT;

					#################
					## SUBROUTINES ##
					#################	

				   #  IN ORDER OF APPEARANCE  #


## READ IN ALIGNMENT FROM FASTA FORMAT 
#
# NOTE SEQ IS RETURNED AS AN ARRAY, SLOW FOR LARGE DATA SETS
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
					#$seq{$name} = $seq;	# STRINGS
					
					@{$seq{$name}} = split (//, $seq); 	# ARRAYS
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

		#$seq{$name} = $seq;				# STRINGS

		@{$seq{$name}} = split (//, $seq); 		# ARRAYS		
	}
	else { print "WARNING: No sequence for $name in fasta file. It is not included in the analysis.\n"; }


	print "\nFound ".@names." sequences in fasta format in $datafile. ";
	if (defined $alignment) { print "Sequences are different lengths\n\n"; }
	else { print "All sequences are $length bp long.\n\n"; }

	return ($length,\%seq,@names);
}


