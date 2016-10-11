#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

my $program = 'fourgamete.pl';					#name of script

my $seq_recog = "0-9A-Z";					# SO INDEXES ARE RECOGNISED
my %parameters;							#input parameters
my $alignment;
my $segfile = "segsites.gfa";
my $maf = 0;
my $secondaf = 0;

getopts('a:eS:m:s:',\%parameters);

if (exists $parameters{"a"}) { $alignment = $parameters{"a"}; }
if (defined $parameters{"S"}) { $segfile = $parameters{"S"}; }
if (defined $parameters{"s"}) { $secondaf = $parameters{"s"}; }
if (defined $parameters{"m"}) { $maf = $parameters{"m"}; }


unless (exists $parameters{"a"}) {
	print "\n USAGE: $program -a <gfa_file>\n\n";
	print   "    -a\tgapped fasta file of aligned seqs\n";
	print   "    -e\texclude sites with > 2 alleles\n";
	print   "    -m\texclude sites with minor allele freq < x [$maf]\n";
	print   "    -s\texclude sites with second allele freq < x [$secondaf]\n";
	print   "    -S\tprint segregating sites to file [$segfile]\n\n";
	die;
}


## READ IN SEQUENCES FROM ALIGNMENT FILE 
#
my ($seq_ref,$length,@names,%seq);
($seq_ref,$length,@names) = fasta2strings($alignment,$seq_recog);	# FASTA FORMAT SEQ
%seq = %$seq_ref;

## IDENTIFY SEGREGATING SITES
#
if (defined $parameters{"s"}) {
	print "\nSites with 2 or more alleles, minor allele frequency > $maf and second allele frequency > $secondaf:\n\n";
}
else {
	print "\nSites with 2 or more alleles and minor allele frequency > $maf:\n\n";
}
my (%segsites, %index);
my $excluded = "no sites excluded";
my $S = 0;
for (my $i=0; $i < $length; $i++) {
	my @bases = qw(A C G T);
	my %alleles;
	my $seqcount = 0;
	foreach my $name (@names) { 
		$seqcount++;
		$seq{$name}[$i] =~ tr/acgtn-/ACGTNN/;
		foreach my $b (@bases) { 
			if ($seq{$name}[$i] eq $b) { $alleles{$b}++; } 
		}
	}

	if (defined $parameters{"s"}) {			# exclude any segregating sites with second alleles with frequencies less than the maf	
		my $allelesovermaf = 0;
		foreach my $a (keys %alleles) {
	#		print "hello2\t$a\t$alleles{$a}\n";
			if ($alleles{$a} > ($secondaf*$seqcount)) { $allelesovermaf++; }
		}
		unless ($allelesovermaf >= 2) { next; }
	}
	my @alleles = keys %alleles;
	if (@alleles < 2) { next; }

	my $minallelecount = $seqcount;
	my $minallele;					# identify the least frequent allele
	foreach my $a (@alleles) {
		if ($alleles{$a} <= $minallelecount) { 	# choose the second one encountered in case of ties
			$minallelecount = $alleles{$a};
			$minallele = $a;
		}
	}
	if ($minallelecount < $maf*$seqcount) { next; }	# excluded segregating sites with < maf

					# include sites with > 2 alleles
	if (!exists $parameters{"e"}) {
		$S++;
		my $pos = $i+1;
		print "$S\t$pos\t".@alleles."\t";	# print summary of sites included
		foreach my $a (@alleles) {
			print "$a: $alleles{$a}\t";
		}
		if (@alleles > 2) { print "$minallele->N\n"; }
		else { print "\n"; }
					# convert the least frequent allele to "N" if there are > 2 alleles at a site
		foreach my $name (@names) { 
			if ((@alleles > 2) && ($seq{$name}[$i] eq $minallele)) {
				push(@{$segsites{$name}},"N");
			}
			else {
				push(@{$segsites{$name}},$seq{$name}[$i]);	
			}
			$index{$S} = $pos;
		}
	}
					# exclude sites with > 2 alleles
	elsif (@alleles == 2) { 
		foreach my $name (@names) { 
			push(@{$segsites{$name}},$seq{$name}[$i]); 
		}
		$S++;
		my $pos = $i+1;
		$index{$S} = $pos;
		print "$S\t$pos\t".@alleles."\t";
		foreach my $a (@alleles) {
			print "$a: $alleles{$a}\t";
		}
		print "\n";
	}
	elsif (@alleles > 2) {
		my $pos = $i+1;
		if ($excluded eq "no sites excluded") { $excluded ="$pos\t".@alleles."\t"; }
		else { $excluded .="$pos\t".@alleles."\t"; }
		foreach my $a (@alleles) {
			$excluded .= "$a: $alleles{$a}\t";
		}
		$excluded .= "\n";
	}

}
if (exists $parameters{"e"}) {
	print "\nWill only consider sites with 2 alleles.\n\n";
	print "Excluded sites:\n$excluded";
}

# PRINT OUT FASTA FILE OF SEGREGATING SITES IF REQUESTED
#
if (exists $parameters{"s"}) {
	open OUT, "> $segfile" or die "couldn't open $segfile : $!";
	foreach my $name (@names) {
		print OUT ">$name\n";
		foreach my $x (@{$segsites{$name}}) { print OUT $x; }
		print OUT "\n";
	}
	close OUT;
	print "\nPrinted a fasta file of segregating sites to $segfile.\n";
}


# IDENTIFY SITES WITH ALL 4 GAMETES
#
print "\nSummary of sites with all 4 gametes:\n";

my (%seen);
my $originalpos = 'none';
for (my $i=0; $i < $S; $i++) {
	for (my $j=0; $j < $S; $j++) {
		if ($i == $j) { next; }

		if (($seen{$i."_$j"}++) || ($seen{$j."_$i"})++) { next; }
		my %gametes;
		foreach my $name (keys %segsites) {
			if ($segsites{$name}[$i] eq "N") { next; }
			if ($segsites{$name}[$j] eq "N") { next; }

			$gametes{"$segsites{$name}[$i]$segsites{$name}[$j]"}++;
		}
		my @diffgametes = keys %gametes;
		unless (@diffgametes == 4) { next; }
		my $ipos = $i+1;
	       	my $jpos = $j+1;	
		print "[$ipos $jpos] ";
		if ($originalpos eq 'none') { $originalpos =  "[$index{$ipos} $index{$jpos}] "; }
		else { $originalpos .=  "[$index{$ipos} $index{$jpos}] "; }


	#	foreach my $g (@diffgametes) {
	#		print "$g: $gametes{$g}\t"
	#	}
	#	print "\n";
 	}

}
print "\n";
print "\nSummary of sites with all 4 gametes in original alignment:\n$originalpos\n\n";



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
					print "Added ".($length-length($seq))." x '-' to $name\n";
					$seq .= '-' x ($length-length($seq));
				}
				if ($seq =~ /[$seq_recog]/i) { 
					unless ($seen{$name}++) { push (@names, $name); }
					else { print "WARNING: $name has been seen $seen{$name} times. Only its last sequence will be analysed.\n"; }
		#			$seq{$name} = $seq;			# STRINGS
					
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

	#	$seq{$name} = $seq;					# STRINGS

		@{$seq{$name}} = split (//, $seq); 			# ARRAYS
	}
	else { print "WARNING: No sequence for $name in fasta file. It is not included in the analysis.\n"; }


	print "\n$datafile :\t".@names." seqs. ";
	if (defined $alignment) { warn "Sequences were different lengths in $datafile\n\n"; }
	else { print "length: $length bp.\n"; }

	return (\%seq,$length,@names);
}


