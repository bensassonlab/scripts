#!/usr/bin/perl


use warnings;
use strict;

open RCMD, ">Rcmds.txt" or die "couldn't open Rcmds.txt : $!";
# load ape
print RCMD "library(ape)\n";
# prepare to construct F84 distance matrix, neighbour joining tree, 
print RCMD "f <- function(x) nj(dist.dna(x,model=\"F84\"))\n";

for (my $i=1; $i<=16; $i++) {
	my $alignment = "CEN$i"."forcbyc.gfa";
	my $outfile = "CEN$i.tree";

	# read alignment
	print RCMD "CEN$i<-read.dna(\"$alignment\",format=\"fasta\")\n";
	# construct F84 distance matrix, neighbour joining tree, 
	# and 10000 (default) bootstraps
	print RCMD "cen$i"."tree<-f(CEN$i)\n";
	print RCMD "set.seed(657)\n"; 		# set seed so that it's repeatable
	print RCMD "cen$i"."bs<-boot.phylo(cen$i"."tree,CEN$i,f,B=10000,quiet=TRUE,trees=TRUE)\n";
	# export tree so you can read it in another program (e.g. figtree)
	print RCMD "cen$i"."tree\$node.label<-cen$i"."bs\$BP\n";
	print RCMD "write.tree(cen$i"."tree,\"$outfile\")\n";
	print RCMD "cen$i"."pp<-prop.part(cen$i"."bs\$trees)\n";
	print RCMD "postprocess.prop.part(cen$i"."pp) # contains bipartitions\n";
	print RCMD "\n";
	print "alignment: $alignment\ttree: $outfile\tbipartitions: Rout.txt\n";
}
close RCMD;

print "\n\nRunning R .. \n";
	
`R CMD BATCH Rcmds.txt Rout.txt`;				# RUN R

print "Done. Phylograms with bootstraps are in newik format in CEN*.tree files\n";


