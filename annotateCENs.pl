#!/usr/bin/perl


use warnings;
use strict;

use Getopt::Std;

my $program = 'annotateCENs.pl';					#name of script

my %parameters;							#input parameters
my $seq_recog = "0-9A-Z\?";					# SO INDEXES ARE RECOGNISED
my ($ref, $file);
my $CDEI = ".TCAC.TG";
my $CDEIII = "G.{8}TTCCGAA.{8}";
my $rCDEI = "CA.GTGA.";
my $rCDEIII = ".{8}TTCGGA[AT].{8}C";
my $outprefix = "CEN";

getopts('i:r:o:',\%parameters);

if (exists $parameters{"i"}) { $file = $parameters{"i"}; }
if (exists $parameters{"r"}) { $ref = $parameters{"r"}; }
if (exists $parameters{"o"}) { $outprefix = $parameters{"o"}; }

unless (exists $parameters{"i"}) {
	print "\n USAGE: $program -i '<seqfile>'\n\n";
	print   "    -i\tseq file in fasta format (gfa, mfa or fa)\n";
	print   "    -r\tname of reference [first sequence]\n";
	print   "    -o\toutput prefix for .gff [$outprefix]\n";
	print "\n";
	die;
}

my ($seq_ref, @names, %seq);
($seq_ref,@names) = fasta2strings($file,$seq_recog);	# FASTA FORMAT SEQ
%seq = %$seq_ref;

unless (defined $ref) {
      	$ref = shift @names;					# assume reference sequence is the first sequence
	print "Reference sequence is first sequence: $ref\n"; 
}
elsif (defined $seq{$ref}) { print "Reference = $ref\n"; }
elsif ((!defined $seq{$ref}) && (defined $parameters{"r"})) {	# or use pattern matching to find it
	foreach my $name (@names) {
		if ($name =~ /$parameters{"r"}/i) { $ref = $name; print "Reference matches $name\n"; }
	}
}
if (!defined $seq{$ref}) { die "ERROR: don't have a sequence for $ref\n"; }


my ($i, $j);
my $cdeIcount = 0;
my $cdeIIIcount = 0;
my $together = 0;

my (@CDEIpos, %together);
while ($seq{$ref} =~ /$CDEI/g) {
	$cdeIcount++;
	$i = pos($seq{$ref});	
	print "CDEI\t".($i-7)."\t$i\t+\n";
       	push(@CDEIpos, $i);	
}
while ($seq{$ref} =~ /$CDEIII/g) {
	$cdeIIIcount++;
	$j = pos($seq{$ref});	
	print "CDEIII\t".($j-23)."\t$j\t+\n";
       	foreach my $cdeIpos (@CDEIpos) {	
		if ($j < ($cdeIpos+200) && ($j > $cdeIpos)) { $together++; $together{$j} = $cdeIpos; }
	}
}
my $s;
if ((defined $i) && (defined $j) && ($together > 10)) { 	# too many hits in +ive strand to search the negative strand
	$s = "+"; 
	foreach my $cdeIIIpos (sort { $a <=> $b } keys %together) {
		print "together\tCEN\t".($together{$cdeIIIpos}-7)."\t$cdeIIIpos\n";
	}
	print "There are so many hits, will cancel searching - strand\n";
}
else {
	if ((defined $i) && (defined $j)) {			# summarise what hits there are then search for more
		$s = "+"; 
		foreach my $cdeIIIpos (sort { $a <=> $b } keys %together) {
			print "together\tCEN\t".($together{$cdeIIIpos}-7)."\t$cdeIIIpos\n";
		}
	}
	undef @CDEIpos;
	if ((defined $i) || (defined $j)) {
		if ($together == 0) { print "WARNING: did not find both CDEI and CDEIII together. Will look in the - strand also.\n"; }
		else { print "Only found $together candidate CENs.  Will look in the - strand also.\n"; }
		undef $i; undef $j; $cdeIcount = 0; $cdeIIIcount = 0;
	}
	while ($seq{$ref} =~ /$rCDEI/g) {
		$cdeIcount++;
		$i = pos($seq{$ref});	
		print "CDEI\t".($i-7)."\t$i\t-\n";
	      	push(@CDEIpos, $i);	
	}
	while ($seq{$ref} =~ /$rCDEIII/g) {
		$cdeIIIcount++;
		$j = pos($seq{$ref});	
		print "CDEIII\t".($j-23)."\t$j\t-\n";
	        foreach my $cdeIpos (@CDEIpos) {	
			if ($j > ($cdeIpos-200) && ($j < $cdeIpos)) { $s = "-"; $together++; $together{$j} = $cdeIpos; }
		}
	
	}
	if ((defined $i) && (defined $j) && (defined $s) && ($s eq "-")) { 
		foreach my $cdeIIIpos (sort { $a <=> $b } keys %together) {
			print "together\tCEN\t".($cdeIIIpos-23)."\t$together{$cdeIIIpos}\n";
		}
	
	
	}
}

if ($together == 0) { die "couldn't find CDEI or CDEIII together in a strand\n"; }
if (($s eq "-") && (($cdeIcount > 1) || ($cdeIIIcount > 1))) { print "WARNING: > 1 hit for CDEI or CDEIII. Might need to add in the correct one by hand\n"; }

open OUT, ">$outprefix.gff" or die "couldn't open $outprefix.gff : $!";
if (($s eq "+") && ($together == 1)) {
	foreach my $cdeIIIpos (keys %together) {
		print OUT $outprefix."CEN\tCEN\t$ref\t".($together{$cdeIIIpos}-7)."\t$cdeIIIpos\t0\t$s\t.\n";
		print OUT $outprefix."CDEI\tCDEI\t$ref\t".($together{$cdeIIIpos}-7)."\t$together{$cdeIIIpos}\t0\t$s\t.\n";
		print OUT $outprefix."CDEII\tCDEII\t$ref\t".($together{$cdeIIIpos}+1)."\t".($cdeIIIpos-24)."\t0\t$s\t.\n";
		print OUT $outprefix."CDEIII\tCDEIII\t$ref\t".($cdeIIIpos-23)."\t$cdeIIIpos\t0\t$s\t.\n";
		print "Printed the 1 CEN annotation (".($together{$cdeIIIpos}-7)."..$cdeIIIpos) to $outprefix.gff\n";
	}
}
elsif ($s eq "+") {
	print OUT $outprefix."CEN\tCEN\t$ref\t".($i-7)."\t$j\t0\t$s\t.\n";
	print OUT $outprefix."CDEI\tCDEI\t$ref\t".($i-7)."\t$i\t0\t$s\t.\n";
	print OUT $outprefix."CDEII\tCDEII\t$ref\t".($i+1)."\t".($j-24)."\t0\t$s\t.\n";
	print OUT $outprefix."CDEIII\tCDEIII\t$ref\t".($j-23)."\t$j\t0\t$s\t.\n";
}
elsif (($s eq "-") && ($together == 1)) {
	foreach my $cdeIIIpos (keys %together) {
		print OUT $outprefix."CEN\tCEN\t$ref\t".($cdeIIIpos-23)."\t$together{$cdeIIIpos}\t0\t$s\t.\n";
		print "Printed the 1 CEN annotation (".($cdeIIIpos-23)."..$together{$cdeIIIpos}) to $outprefix.gff\n";
	}
}
elsif ($s eq "-") {
	print OUT $outprefix."CEN\tCEN\t$ref\t".($j-23)."\t$i\t0\t$s\t.\n";
	print OUT $outprefix."CDEI\tCDEI\t$ref\t".($i-7)."\t$i\t0\t$s\t.\n";
	print OUT $outprefix."CDEII\tCDEII\t$ref\t".($j+1)."\t".($i-8)."\t0\t$s\t.\n";
	print OUT $outprefix."CDEIII\tCDEIII\t$ref\t".($j-23)."\t$j\t0\t$s\t.\n";
}
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


