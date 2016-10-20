#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

my $program = 'fa2phylip.pl';					#name of script

my $seq_recog = "0-9A-Z";					# SO INDEXES ARE RECOGNISED
my $block_length = 100;
my (%parameters, $in);						#input parameters

getopts('i:qb:T',\%parameters);


if (exists($parameters{"i"})) { $in = $parameters{"i"}; }
if (exists($parameters{"b"})) { $block_length = $parameters{"b"}; }
if (exists($parameters{"T"})) { $block_length = 10000; }


unless (exists $parameters{"i"}) {
	print "\n USAGE: $program -i <fasta_file>\n\n";
	print " -q\tRun in quiet mode: no warnings if sequences are excluded.\n";
	print " -b\tBlocklength [$block_length].\n";
	print " -T\tFor use by TCS\n\n";
	print " All sequences must be of equal length.\n\n";
	exit;
}


## READ IN SEQUENCES FROM ALIGNMENT FILE 
#
my ($seq_ref,@names,%seq);
($seq_ref,@names) = fasta2strings($in,$seq_recog);	# FASTA FORMAT SEQ
%seq = %$seq_ref;

## MAKE NAMES THE SAME LENGTH BY ADDING BLANKS OR TRUNCATING NAME
my $max = 10; my %newname;
foreach my $name (@names) { if ($name =~ /(\S{$max})\S+/) { $newname{$name} = $1; }}	# truncate name to <=10 char	
foreach my $name (@names) { 
	if (exists $parameters{"T"}) {
		my $newname;
		if (length($name) > 9) { ($newname = $name) =~ s/_//g; $newname =~ /^(.{9})/m; $newname{$name} = $1; }
		if (length($name) < 9) { $newname{$name} = $name. "x" x (9-length($name)); } 
	}
	elsif (length($name) < $max) { 
		if (exists $parameters{"T"}) { $newname{$name} = $name. "x" x ($max-length($name)); } 
		else { $newname{$name} = $name. " " x ($max-length($name)); } 
	}
}



## PRINT OUT ALIGNMENT IN PHYLIP FORMAT
#
print @names."  ".length($seq{$names[0]})."\n";
foreach my $name (@names) {
	my $seq = substr($seq{$name},0,$block_length);
	if (exists $newname{$name}) { print "$newname{$name} $seq\n"; }
	else { print "$name $seq\n"; }

}
print "\n";

for (my $i=$block_length; $i <= length($seq{$names[0]}); $i+=$block_length) {
	foreach my $name (@names) {
		my $seq = substr($seq{$name},$i,$block_length);
		print "$seq\n";
	}
	print "\n";
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
				else { unless (exists $parameters{'q'}) { print "WARNING: No sequence for $name in fasta file. It is not included in the analysis.\n"; } }
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
	else { unless (exists $parameters{'q'}) { print "WARNING: No sequence for $name in fasta file. It is not included in the analysis.\n"; } }

	if (defined $alignment) { die "Sequences are different lengths\n\n"; }

	return (\%seq,@names);
}

