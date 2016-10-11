#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

# note: I tested this with the following non-default structure settings in mainparams
# #define PLOIDY       1    // (int) ploidy of data
# #define ONEROWPERIND 1    // (B) store data for individuals in a single line
# #define MARKERNAMES      0  // (B) data file contains row of marker names
# #define MAPDISTANCES     1  // (B) data file contains row of map distances

my $program = 'structureShell.pl';				#name of script

my %parameters;							#input parameters
my $infile;
my $L = 53;
my $N = 83;
my $maxK = 15;
my $minK = 11;

getopts('i:L:N:s:e:',\%parameters);

if (exists $parameters{"i"}) { $infile = $parameters{"i"}; }
if (exists $parameters{"L"}) { $L = $parameters{"L"}; }
if (exists $parameters{"N"}) { $N = $parameters{"N"}; }
if (exists $parameters{"s"}) { $minK = $parameters{"s"}; }
if (exists $parameters{"e"}) { $maxK = $parameters{"e"}; }

unless (exists $parameters{"i"}) {
	print "\n USAGE: $program -i <infile>\n\n";

	print   "    -i\tinfile for structure\n";
	print   "    -s\tstarting K [$minK]\n";
	print   "    -e\tending K [$maxK]\n";
	print   "    -L\tnumber of loci [$L]\n";
	print   "    -N\tnumber of individuals [$N]\n\n";
	exit;
}

for (my $k=$minK; $k<=$maxK; $k++) {
	open OUT, ">$infile.k$k.structure.stdout" or die "couldn't open $infile.k$k.structure : $!";
	print "\nK = $k\n"; 
	open MAINPARAMS, "<mainparams" or die "couldn' open mainparams : $!";
	open TEMP, ">temp" or die "couldn't open temp : $!"; 
	while (<MAINPARAMS>) {
		if (/#define MAXPOPS/) { print TEMP "#define MAXPOPS    $k      // (int) number of populations assumed\n"; }
		if (/#define MAXPOPS/) { print "#define MAXPOPS    $k      // (int) number of populations assumed\n"; }
		else { print TEMP $_; }
	}
	close TEMP;
	close MAINPARAMS;
#	system ("cat mainparams");
	`cp temp mainparams`;

	my $structureout = `structure -i $infile -L $L -N $N`;

	while ($structureout =~ /Estimated Ln Prob of Data   = (\S+)/g) {
		print "Estimated Ln Prob of Data   = $1\n";
	}
 	print OUT $structureout;
	close OUT;
	`cp outfile_f $infile.k$k.structure`;
}
