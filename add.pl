#!/usr/bin/perl

use strict;

my $program = 'add.pl';					# name of script
my $total;

open (IN, "-");
while (<IN>) {
	/([\d\.-]+)/;
	$total += $1; 
}
print "$total\n";


