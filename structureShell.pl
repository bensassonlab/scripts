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
my $replicates = 1;

getopts('i:L:N:s:e:r:S',\%parameters);

if (exists $parameters{"i"}) { $infile = $parameters{"i"}; }
if (exists $parameters{"L"}) { $L = $parameters{"L"}; }
if (exists $parameters{"N"}) { $N = $parameters{"N"}; }
if (exists $parameters{"s"}) { $minK = $parameters{"s"}; }
if (exists $parameters{"e"}) { $maxK = $parameters{"e"}; }
if (exists $parameters{"r"}) { $replicates = $parameters{"r"}; }



unless (exists $parameters{"i"}) {
	print "\n USAGE: $program -i <infile>\n\n";

	print   "    -i\tinfile for structure\n";
	print   "    -s\tstarting K [$minK]\n";
	print   "    -e\tending K [$maxK]\n";
	print   "    -r\tnumber of times to repeat STRUCTURE run from starting K to ending K [$replicates]\n";
	print   "    -S\tuse with -r to change seed for each replicate run of STRUCTURE\n";
	print   "    -L\tnumber of loci [$L]\n";
	print   "    -N\tnumber of individuals [$N]\n\n";
	exit;
}

for (my $r=1; $r<=$replicates; $r++) {			# run structure $replicates times 

	my $seed = int(rand(10000));

	for (my $k=$minK; $k<=$maxK; $k++) {		# from $minK to $maxK
		my $outfile;
		if ($replicates > 1) { 
			if (defined $parameters{'s'}) { $outfile = "$infile.k$k"."_$r"."_$seed.structure"; }
			else { $outfile = "$infile.k$k"."_$r.structure"; }
		}
		else { $outfile = "$infile.k$k.structure"; }
		open OUT, ">$outfile.stdout" or die "couldn't open $outfile.stdout : $!";
		print "\nK = $k\n";
		if ($replicates > 1) { print "Replicate $r\n"; }
		open MAINPARAMS, "<mainparams" or die "couldn' open mainparams : $!";	# change K (MAXPOPS) in mainparams
		open TEMP, ">temp" or die "couldn't open temp : $!"; 
		while (<MAINPARAMS>) {
			if (/#define MAXPOPS/) { print TEMP "#define MAXPOPS    $k      // (int) number of populations assumed\n"; }
			if (/#define MAXPOPS/) { print "#define MAXPOPS    $k      // (int) number of populations assumed\n"; }
			else { print TEMP $_; }
		}
		close TEMP;
		close MAINPARAMS;

		`cp temp mainparams`;

		if (defined $parameters{'s'}) {
			open EXTRAPARAMS, "<extraparams" or die "couldn' open extraparams : $!";	# change seed in extraparams
			open TEMP, ">temp" or die "couldn't open temp : $!"; 
		
			while (<EXTRAPARAMS>) {

				if (/#define SEED/) { print TEMP "#define SEED    $seed      //  (int) seed value for random number generator perl: int(rand(10000))\n"; }
				if (/#define SEED/) { print "#define SEED    $seed      // (int) seed value for random number generator perl: int(rand(10000))\n"; }
				else { print TEMP $_; }
			}
			close TEMP;
			close EXTRAPARAMS;

			`cp temp extraparams`;
		}

		my $structureout = `structure -i $infile -L $L -N $N`;

		while ($structureout =~ /Estimated Ln Prob of Data   = (\S+)/g) {
			print "Estimated Ln Prob of Data   = $1\n";
		}
	 	print OUT $structureout;
		close OUT;
		`cp outfile_f $outfile`;
	}
}
