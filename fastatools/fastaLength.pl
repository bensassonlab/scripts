#! /usr/bin/perl

use strict;
use Getopt::Std;						# read in data at commandline using -i type options

my $program = 'fastaLength.pl';					# name of script

my (%parameters, @names, $name, %seq);


### INPUT ###

getopts('i:gN',\%parameters);

unless (exists $parameters{"i"}) {
	print "\n\n\nUSAGE: $program -i <fasta_seqs>\n\n";
	print "\t-g  show ungapped length\n";
	print "\t-N  count Ns as gaps\n\n";
}


open IN, $parameters{'i'} or die "couldn't open $parameters{'i'} : $!";

while (<IN>) {						# read fasta %seq
	if (/>(.*)/) { $name = $1; push (@names, $name); }	# retain alignment order
	elsif (/^([a-z~\?-]+)\s*$/mi) { $seq{$name} .= $1; } 
}
print "\nFound ".@names." sequences in $parameters{'i'}.\n\n";


foreach (@names) { 
	print "$_\t".length($seq{$_})." bp";
	$seq{$_} =~ tr/~\?-//d;
	if (exists $parameters{"N"}) { $seq{$_} =~ tr/Nn//d; }
	unless (exists $parameters{'g'}) { 
		unless (exists $parameters{'N'}) { print "\n"; }
		else {	print " ungapped[Nn] ".length($seq{$_})." bp\n"; }
	}
	else {
	       if (exists $parameters{'N'}) {
		       	print " ungapped[~Nn-?] ".length($seq{$_})." bp\n"; 
	       }
	       else {	print " ungapped[~-?] ".length($seq{$_})." bp\n"; }
	} 
}
