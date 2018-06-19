#!/usr/bin/perl

# this is a test line to see how github repo updates on kadmon

use warnings;
use strict;
use Getopt::Std;

my $program = 'basecomp.pl';					#name of script

my $seq_recog = "acgtnxykwrsm";					# MOST AMBIGUITY CODES ARE RECOGNISED FOR MOST SEQUENCES
my @bases = ('A','C','G','T');
my %exclude;
$exclude{'-'} = 1;

my %parameters;							#input file names
my ($datafile, @names, %seq, $bases, %count, %total_bases, %sum, %allcharacters);		
my $total = 0;


getopts('i:n:sp',\%parameters);


if (exists($parameters{"i"})) { $datafile = $parameters{"i"}; }
if (exists($parameters{"n"})) { $datafile = $parameters{"n"}; }


unless ((exists $parameters{"i"}) || (exists $parameters{"n"})) {
	print "\n USAGE: $program -i <datafile>\n";
	print   " OR:    $program -n <datafile>\n\n"; 
	print   "    -i\tdatafile: sequence alignment in fasta format\n";
	print 	"    -n\tsequence alignment in nexus format\n";
	print 	"    -p\tshow proportions not counts\n";
	print 	"    -s\tsummarise with overall totals or mean proportions\n\n";
	exit;
}


## READ IN SEQUENCES FROM FILE 
#
my $seq_ref;

if (exists $parameters{'n'}) { ($seq_ref,@names) = nexus($datafile); }	# NEXUS FORMAT
else { ($seq_ref,@names) = fasta2strings($datafile,$seq_recog); }	# FASTA FORMAT

%seq = %$seq_ref;


foreach my $base (@bases) { $bases .= $base; }				# LOOKUP STRING FOR BASES

foreach my $name (@names) {

	while ($seq{$name} =~ /(\S)/gi) {				# GET COUNT INFO FOR ALL CHARACTERS
		my $character = $1;
		$allcharacters{$character}++;
		$count{"$name$character"}++;
	}		
	
	foreach my $base (@bases) { 
#		while ($seq{$name} =~ /$base/gi) {			# COUNT EVERY SPECIAL BASE: CASE INSENSITIVE
#			$count{"$name$base"}++;
#		}	
		if (exists $parameters{'p'}) { 
			$total_bases{$name} += $count{"$name$base"}; 
			$total_bases{'all'} += $count{"$name$base"}; 
		}
		if (exists $parameters{'s'}) { $sum{$base} += $count{"$name$base"}; }
		
	}
}	

## PRINT SUMMARY (total all bases)
#
if (exists $parameters{'s'}) {
	foreach my $base (@bases) { 					# SUMMARY OF SPECIAL BASES
		$total += $sum{$base};
		if (exists $parameters{'p'}) { print "\t$base : ".($sum{$base}/$total_bases{'all'}); }
		else { print "\t$base : $sum{$base}"; }
	}
	print "\tTotal$bases: $total";

       foreach my $other (sort keys %allcharacters) {			# PRINT OUT A SUMMARY OF OTHER CHARACTERS
	       	if (($bases !~ /$other/) && (!defined $exclude{$other})) {
			$total += $allcharacters{$other};
			print "\t$other : $allcharacters{$other}";
		}
	}
	print "\tTotalALL: $total";
		
}
## PRINT RESULTS FOR EVERY SEQUENCE
#
else {
	foreach my $name (@names) {
		my $ntotal = 0;
		print "$name"; 
	       	foreach my $base (@bases) { 				# SUMMARY OF SPECIAL BASES
			if(defined $count{"$name$base"}) { $ntotal += $count{"$name$base"}; }
			unless (exists $count{"$name$base"}) { $count{"$name$base"} = 0; }
			if (exists $parameters{'p'}) {  
				print "\t$base : ".($count{"$name$base"}/$total_bases{$name});
			}
			else { print "\t$base : ".$count{"$name$base"}; }
	       }
	       print "\ttotal$bases : $ntotal";
	       
	       foreach my $other (sort keys %allcharacters) {		# PRINT OUT A SUMMARY OF OTHER CHARACTERS
		       if (($bases !~ /$other/) && (!defined $exclude{$other})) {
		       		unless (exists $count{"$name$other"}) { $count{"$name$other"} = 0; }
				$ntotal += $count{"$name$other"};
		       		print "\t$other : ".$count{"$name$other"};
		       }
	       }
	       print "\ttotalALL : $ntotal";

	       unless (length($seq{$name}) == $ntotal) { 
		       die "ERROR: the total number of bases recognised does not equal the length of sequence."; 
	       }

       	       print "\n";	       
	} 	
}
print "\n";


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

	return (\%seq,@names);
}

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
		elsif (/^\s*([$seq_recog\-]+)\s*$/mi) { $seq .= $1; }	# GATHER SEQ (even if spread across multiple lines)
		elsif (/\S+/) { 
			while (/(\S)/g) { 
				unless ($1 =~ /([$seq_recog\-])/i) { 
					print "seq_recog[$seq_recog] unrecognised character: [$1]\n"; 
				} 
			}
		}
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

