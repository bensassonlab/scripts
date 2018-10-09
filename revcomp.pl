#!/usr/bin/perl

use strict;

my $program = 'revcomp.pl';					# name of script
my $revcomp_seq;

open (IN, "-");
while (<IN>) {
	if (/^>/m) { print $_; next; } 
	/([a-z]+)/i;
	$revcomp_seq = reverse $1;
	$revcomp_seq =~ tr/acgtACGT/tgcaTGCA/;
	print "$revcomp_seq\n";
}


