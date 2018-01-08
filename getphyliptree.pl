#!/usr/bin/perl


use warnings;
use strict;

use Getopt::Std;

my $program = 'getphyliptree.pl';					#name of script

my %parameters;							#input parameters
my $gfafile;
my $n = 10;

getopts('i:n:',\%parameters);

if (exists $parameters{"i"}) { $gfafile = $parameters{"i"}; }
if (exists $parameters{"n"}) { $n = $parameters{"n"}; }


unless (exists $parameters{"i"}) {
	print "\n USAGE: $program -i '<fastafile>'\n\n";
	print   "    -i\tfasta format alignment\n";
	print   "    -n\tnumber of bootstrap replicates [$n]\n\n";
	print   " NOTE: Requires fa2phylip.pl in your path\n\n";
	exit;
}

								# CONVERTING GFA TO PHY
`fa2phylip.pl -i $gfafile > $gfafile.phy`;

print "\nRunning neighbour joining distance analysis with $n bootstraps to generate a consensus tree ..\n\n";

								# RUNNING DNAML
#open DNAMLIN, ">dnamlin" or die "couldn't open dnamlin: $!";
#print DNAMLIN "$gfafile.phy\nT\n0.5\nY\n";
#`dnaml < dnamlin > dnamlout`;
#`mv outfile $gfafile.dnamlout`;
#`mv outtree $gfafile.mltree.phy`;

								# RUNNING DNADIST
open DNADISTIN, ">dnadistin" or die "couldn't open dnadistin: $!";
print DNADISTIN "$gfafile.phy\nT\n0.5\nY\n";
`dnadist < dnadistin > dnadistout`;
`mv outfile $gfafile.dnadistout`;

								# RUNNING NEIGHBOR
open NEIGHBORIN,  ">neighborin" or die "couldn't open neighborin: $!";
print NEIGHBORIN "$gfafile.dnadistout\nJ\n3\nY\n";
`neighbor < neighborin > neighborout`;
`mv outfile $gfafile.neighborout`;
`mv outtree $gfafile.tree.phy`;

								# RUNNING SEQBOOT
open SEQBOOTIN,  ">seqbootin" or die "couldn't open seqbootin: $!";
print SEQBOOTIN "$gfafile.phy\nR\n$n\nY\n3\n";
`seqboot < seqbootin > seqbootout`;
`mv outfile $gfafile.seqbootout`;

								# RUNNING DNAML FOR BOOTSTRAPS
#open DNAMLIN, ">dnamlin" or die "couldn't open dnamlin: $!";
#print DNAMLIN "$gfafile.seqbootout\nM\nD\n$n\n3\n5\nT\n0.5\nY\n";
#`dnaml < dnamlin > dnamlBootout`;
#`mv outfile $gfafile.dnamlBootout`;
#`mv outtree $gfafile.mltrees.phy`;


								# RUNNING DNADIST FOR BOOTSTRAPS
open DNADISTIN, ">dnadistin" or die "couldn't open dnadistin: $!";
print DNADISTIN "$gfafile.seqbootout\nT\n0.5\nM\nD\n100\nY\n";
`dnadist < dnadistin > dnadistBootout`;
`mv outfile $gfafile.dnadistBootout`;

								# RUNNING NEIGHBOR FOR BOOTSTRAPS
open NEIGHBORIN,  ">neighborin" or die "couldn't open neighborin: $!";
print NEIGHBORIN "$gfafile.dnadistBootout\nJ\n3\nM\n100\n3\nY\n";
`neighbor < neighborin > neighborBootout`;
`mv outfile $gfafile.neighborBootout`;
`mv outtree $gfafile.trees.phy`;

								# RUNNING CONSENSE
open CONSENSEIN,  ">consensein" or die "couldn't open consensein: $!";
# print CONSENSEIN "$gfafile.mltrees.phy\nY\n";			# for MAX LIKE
print CONSENSEIN "$gfafile.trees.phy\nY\n";			# for neighbour joining
`consense < consensein > consenseout`;
`mv outfile $gfafile.consenseout`;
`mv outtree $gfafile.cons.tree`;

print "Neighbour joining tree is in $gfafile.tree.phy.\n";

print "Final consensus tree of the bootstrap data is in $gfafile.cons.tree, \nand summary of bootstrap support is in $gfafile.consenseout\n\n";

