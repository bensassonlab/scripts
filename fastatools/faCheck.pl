#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

my $program = 'faCheck.pl';					#name of script

my $seq_recog = "0-9A-Z\?";					# SO INDEXES ARE RECOGNISED
my %parameters;							#input parameters
my ($seqfile1, $seqfile2);

getopts('i:j:na',\%parameters);


if (exists($parameters{"i"})) { $seqfile1 = $parameters{"i"}; }
if (exists($parameters{"j"})) { $seqfile2 = $parameters{"j"}; }


unless ((exists $parameters{"i"}) && (exists $parameters{"j"})) {
	print "\n USAGE: $program -i <fasta_seqs> -j <fasta_seqs>\n\n";
	print   "    -i\tfasta file of seqs\n";
	print 	"    -j\tfasta file of seqs to compare to i\n";
	print 	"    -n\tdo not treat Ns as differences\n";
	print 	"    -a\tignore N insertions: this needs alignments to be same length\n\n";
	print 	" Note: all 'X' and '-' differences are ignored\n\n";
	exit;
}


## READ IN SEQUENCES FROM ALIGNMENT FILE 
#
my ($seq_ref,@names1,@names2,%seq1,%seq2);
($seq_ref,@names1) = fasta2strings($seqfile1,$seq_recog);	# FASTA FORMAT SEQ
%seq1 = %$seq_ref;

($seq_ref,@names2) = fasta2strings($seqfile2,$seq_recog);	# FASTA FORMAT SEQ
%seq2 = %$seq_ref;

for (my $i=0; $i < @names1; $i++) {
	unless (exists $seq2{$names1[$i]}) { 
		print "$seqfile1\t$names1[$i]\tmissing in $seqfile2\n";
		next;	
	}
	my ($seq1,$seq2);
	$seq1{$names1[$i]} =~ tr/[a-z]/[A-Z]/;
	$seq2{$names1[$i]} =~ tr/[a-z]/[A-Z]/;
	($seq1 = $seq1{$names1[$i]}) =~ tr/X-//d;		
	($seq2 = $seq2{$names1[$i]}) =~ tr/X-//d; 
	if ($names1[$i] eq 'spar_chr3') {
		$seq1 =~ tr/Nn//d;
		$seq2 =~ tr/Nn//d; 
	}		
	if (exists $parameters{'n'}) {
		$seq1 =~ tr/Nn/../;
		$seq2 =~ tr/Nn/../;
		#	die "\n\n$seq1\n";
	}

	my $compared = 0;
	if ($seq1 eq $seq2) {
	       	if (length($seq2) == length($seq1)) { 
			$compared = length($seq2); 	
			print "$seqfile1\t$names1[$i]\tidentical in $seqfile2\tcompared: $compared bp\n";
		}	
		else { print "Error: seqs are identical, but different lengths ".length($seq1)." and ".length($seq2)."\n"; }
	}
	else { 
		unless (exists $parameters{'n'}) {
			print "$seqfile1\t$names1[$i]\t".length($seq1)." bp\tdifferent in $seqfile2\t".length($seq2)." bp\n";
		}

		my @seq1 = split (//,$seq1{$names1[$i]});
		my @seq2 = split (//,$seq2{$names1[$i]});
		my @gappedseq1 = @seq1;			# ???
		my @gappedseq2 = @seq2;			# ???

		my ($diff1,$diff2,$j,%seq2gapped1,%seq2gapped2);
		my $count = 0;
		for (my $i = 0; $i < @seq1; $i++) {
			if (!defined $j) { $j = 0; }
			else { $j++; }
			if (!defined $seq2[$j]) { last; }


			while ($seq1[$i] =~ /[X-]/g) { $i++; if ($i == @seq1) { last; } }
			while ($seq2[$j] =~ /[X-]/g) { $j++; if ($j == @seq2) { last; } } 
			if ($i == @seq1) { last; } 	
			
#			if ($i != $j) { print "$seqfile1 $i $seq1[$i] ne $seqfile2 $j $seq2[$j]\n"; }		
			
			if (($seq1[$i] ne $seq2[$j]) || (($count > 0) && ($count < 20))) {
				if ((exists $parameters{'n'}) && (($seq1[$i] =~/n/i) || ($seq2[$j] =~ /n/i))) { next; }	# ignore differences due to "n"s

				if (exists $parameters{'a'}) {

				# ignore differences if bases in the original alignments were the same (this way insertions caused by Ns 
				# can still be ignored if alignment coordinates are exactly equivalent between the 2 alignments)

					if ((exists $parameters{'n'}) && ($gappedseq1[$i] eq $gappedseq2[$i])) { next; }
					if ((exists $parameters{'n'}) && ($gappedseq1[$j] eq $gappedseq2[$j])) { next; }		# ???
					if ((exists $parameters{'n'}) && ($gappedseq1[$i] =~ /n/i)) { next; }
					if ((exists $parameters{'n'}) && ($gappedseq2[$i] =~ /n/i)) { next; }		# ???
					if ((exists $parameters{'n'}) && ($gappedseq1[$j] =~ /n/i)) { next; }
					if ((exists $parameters{'n'}) && ($gappedseq2[$j] =~ /n/i)) { next; }		# ???

				}

				if (($count == 0) && (exists $parameters{'n'})) {
					print "$seqfile1\t$names1[$i]\t".length($seq1)." bp\tdifferent in $seqfile2\t".length($seq2)." bp (wild Ns)\n";
				}

				print "$seq1[$i] pos ".($i+1)." in $seqfile1 ne $seq2[$j] pos ".($j+1)." in $seqfile2\n";
				$diff1 .= $seq1[$i];
			       	$diff2 .= $seq2[$j];
				$count++;
			}
			if ($count >= 20) {
			       	print "Diffs start with : \t$diff1 in $seqfile1\n";
				print "                   \t$diff2 in $seqfile2\n";
				last; 
			}
			unless ( ($gappedseq1[$i] =~ /[nX-]/i) || ($gappedseq2[$i] =~ /[nX-]/i) ) { $compared++; }	
		}
		if ((exists $parameters{'n'}) && ($count == 0)) { 
			print "$seqfile1 $names1[$i]\tsame in $seqfile2 ".length($seq2)." bp (wild Ns) compared: $compared\n"; 		
		}
		elsif ((exists $parameters{'n'}) && ($count > 0)) { 
			print "$seqfile1 $names1[$i]\tdifferent in $seqfile2 ".length($seq2)." bp (wild Ns) WARNING compared: over$compared\n"; 
		}


	}      	
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

