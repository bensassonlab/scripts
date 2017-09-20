#!/usr/bin/perl

use strict;

my $program = 'subreadLcheck.pl';					# name of script
my $totalNo = 0;
my $over10kb = 0;
my $over15kb = 0;
my $over20kb = 0;
my $over30kb = 0;

my $totalLength = 0;
my $over10kblength = 0;
my $over15kblength = 0;
my $over20kblength = 0;
my $over30kblength = 0;


print "\n# READ COUNTS\ntotalNo | over10kb | over15kb | over20kb | over30kb\n";
print " -------| ---------| -------- | ---------| --------\n";

open (IN, "-");
while (<IN>) {
	/([\d\.-]+)/;
	my $length = $1;
	$totalNo++; $totalLength += $length;
	if ($length > 10000) { $over10kb++; $over10kblength+=$length; }
	if ($length > 15000) { $over15kb++; $over15kblength+=$length; }
	if ($length > 20000) { $over20kb++; $over20kblength+=$length; }
	if ($length > 30000) { $over30kb++; $over30kblength+=$length; }

}
print "$totalNo | $over10kb | $over15kb | $over20kb | $over30kb\n";

print "\n# TOTAL LENGTH (Gbp)\ntotalLength | over10kb | over15kb | over20kb | over30kb\n";
print " -------| ---------| -------- | ---------| --------\n";

$totalLength /= 1000000000;
$over10kblength /= 1000000000;
$over15kblength /= 1000000000;
$over20kblength /= 1000000000;
$over30kblength /= 1000000000;


print "$totalLength | $over10kblength | $over15kblength | $over20kblength | $over30kblength\n";



