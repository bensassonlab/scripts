#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

# note: I tested this with the following non-default structure settings in mainparams
# #define PLOIDY       1    // (int) ploidy of data
# #define ONEROWPERIND 1    // (B) store data for individuals in a single line
# #define MARKERNAMES      0  // (B) data file contains row of marker names
# #define MAPDISTANCES     1  // (B) data file contains row of map distances

my $program = 'structurePrint.pl';				#name of script

my %parameters;							#input parameters
my ($infile,@un,$datafile);
my $NoOfcolours = 10;
my $seed = 54;
my $start = 0;
my $end = 1;
my $colours = '"#8c510a","#4daf4a","#377eb8","#e41a1c","orange","#ffff33"';


print "\nThis script uses R to create a plot from a STRUCTURE outfile.\n\nDependencies: R\n\n";

getopts('i:n:c:C:s:e:S:lq',\%parameters);

if (exists $parameters{"i"}) { $infile = $parameters{"i"}; $datafile = "$infile.data"; }
if (defined $parameters{"c"}) { $NoOfcolours = $parameters{"c"}; }
if (defined $parameters{"C"}) { $colours = $parameters{"C"}; }
if (defined $parameters{"s"}) { $start = $parameters{"s"}; }
if (defined $parameters{"e"}) { $end = $parameters{"e"}; }
if (defined $parameters{"S"}) { $seed = $parameters{"S"}; }
if (exists $parameters{"n"}) { 			# user names provided by user to define popn print order
	@un = split(/\s+/, $parameters{"n"});

# note structure truncates names to 11 characters e.g. UWOPS03_461 
# therefore need to truncate user names that way too
	my @temp;
	foreach my $name (@un) {
		if ($name =~ /^(\S{11})\S+/m) { push(@temp, $1); }
		else { push (@temp, $name); }
	}
	@un = @temp; 	
	print "names that define populations: @un\n\n";	
}

unless (exists $parameters{"i"}) {
	print "\n USAGE: $program -i <infile>\n\n";

	print   "    -i\tstructure outfile e.g. outfile_f\n";
	print   "    -q\tsort by Q (similar to STRUCTURE)\n";
	print   "    -n\tnames that define popn order e.g. 'YJM981m YPS606m'\n";
	print   "    -l\tdraw with name labels shown on columns\n";
	print   "    -C\tcolour choices e.g. '$colours'\n";
	print   "    -c\tmaximum colours to use for plots [$NoOfcolours]\n";
	print   "    -s\tstart of colours - see rainbow in R [$start]\n";
	print   "    -e\tend of colours - see rainbow in R [$end]\n";
	print   "    -S\tseed for colours - see rainbow in R [$seed]\n\n";


	exit;
}
open IN, "<$infile" or die "couldn't open $infile : $!";
open DATA, ">$datafile" or die "couldn't open $datafile : $!"; 

my $k = 0;
my $pk;				# find out K for printing from USEPOPINFO structure outfiles

my (%data,@names,%rows,@pops,%seenPop,%maxq);
while (<IN>) {

       	if (/(\d+) populations assumed/) { $pk = $1; }	
		
	if (/^(\s*\d+\s+(\S+)\s+\(\d+\)\s+(\d+)\s*:\s*)(\S+.*)\s*$/m) {
		my $start = $1;
		my $name = $2;
		my $pop = $3; 		# find out user-defined pop for printing from USEPOPINFO structure outfiles
		my $p = $4;
		$k = 0;

		if ($p =~ /\*\*\*/) { 	# assume user defined USEPOPINFO individuals are solely from user-specified population
			my $newp;
			for (my $i=0; $i<$pk; $i++) {
				if ($i+1 == $pop) { $newp .= "1.000\t"; }
				else { $newp .= "0.000\t"; }
			}
			$p = $newp;
			$name .= "***";	# mark name, so that it is clear this is a user-specified individual
			print "$name\tuser-specified\tpopn: $pop\n";
		}	

		push (@names,$name);

		while ($p =~ /(\S+)/g) { 
			$k++; 
			push(@{$data{$k}},$1);	# gather up data by popNo to decide sort order
			push(@{$data{$name}},$1);	# gather up data by Individual to decide sort order
			unless($seenPop{$k}++) { push(@pops, $k); }
		} 
		my @line = split (/\s+/,"$start$p"); 
		@{$rows{$name}} = @line;
						# figure out the maximum value of q (ie the best popn assignment)
		for (my $qi=6; $qi < @line; $qi++) {
			if (!defined $maxq{$name}) { $maxq{$name} = $line[$qi]; }
			elsif ($line[$qi] > $maxq{$name}) { $maxq{$name} = $line[$qi]; }
		}
		unless (exists $parameters{"q"}) { 
			foreach my $d (@line) { print DATA "$d\t"; } 
			print DATA "\n";	
		}
	}
}

if (($k > $NoOfcolours) && (!defined $parameters{"c"})) {
       print "warning: > 1 population may have the same colour if k>12. Please set the number of colours using the-c option.\n"; 
}

my (@csorted,%printed);
my $check = 0;
if (exists $parameters{"q"}) {			# difficult sort, but looks good. 
						# Easier for me to implement in perl than in R
	my (@poporder,%seen);
	foreach my $n (@un) {
		if (defined @{$data{"$n***"}}[0]) { $n .= "***"; }
		if (!defined @{$data{$n}}[0]) { print "ERROR: $n is not in structure outfile: $infile\n"; next; }
		my $max = 0;
		my $pop;			# which popns do the user-named individuals belong to?
		for (my $i=0; $i < @{$data{$n}}; $i++) { 
			if ($data{$n}[$i] > $max) { $max = $data{$n}[$i]; $pop = $i+1; }
		}
		unless ($seen{$pop}++) { push(@poporder,$pop); }
		print "$n\tpop=$pop\t@{$data{$n}}\n";

	}
	foreach my $pop (@pops) {
		if (!defined $seen{$pop}) { push (@poporder,$pop); }
	}

	foreach my $c (@poporder) { 	# print q values in order of k columns in structure results file
						# sort popns for print using names from user	
		my @q =  @{$data{$c}}; 
		my @sid = sort { $q[$b] <=> $q[$a] } 0 .. $#q;	# an index \@sid for sorting multiple arrays

		my @sq =  @q[@sid];
		my @snames = @names[@sid];

	
		for (my $i=0; $i < @sq; $i++) {
							# when q is below 0.5 if this is not the best popn
							# assignment, then move onto next cluster

			if (($sq[$i] < 0.5) && ($sq[$i] < $maxq{$snames[$i]})) { 
			       	next; 
			}									
			$check++;
			for (my $j=1; $j<=5; $j++) {
				if ($j == 2) { print DATA "$snames[$i]\t"; next; }
				print DATA "@{$rows{$snames[$i]}}[$j]\t";
			}
		#	my @sortedq = sort { $a <=> $b } @{$data{$snames[$i]}};
		#	print DATA "@sortedq\n";
			foreach my $j (@poporder) {
				print DATA @{$rows{$snames[$i]}}[$j+5]."\t"; 
			}
			print DATA "\n"
		}
	
		
	}

}
if ($check != @names) { print "WARNING: $check lines printed to $datafile not ".@names."\n"; }
close DATA;


print "\n\n$infile\tK = $k\n";

open RCMD, ">Rcmds.txt" or die "couldn't open Rcmds.txt : $!";
print RCMD "data<-read.table(\"$datafile\")\n";
print RCMD "attach(data)\n";

print RCMD "x<-rbind(data[,6]"; 			# create a matrix of data columns showing proportions
for (my $i=7; $i<(6+$k); $i++) { print RCMD ",data[,$i]"; }
print RCMD ")\n";

print RCMD "pdf('$infile.pdf',width=10, height=4)\n";	# name plotfile and set the shape of the plot

#print RCMD "par(mfrow=c(3,1))\n";
#for (my $seed=1; $seed<=100; $seed++) {
print RCMD "par(mar=c(7, 4, 4, 2) + 0.1,cex=0.8)\n";
if (defined $parameters{"C"}) { 			# User specified colours
	print RCMD "barplot(x,beside=FALSE,col=c($colours),"; 
} 
else {							# random colors
	print RCMD "set.seed($seed)\n";			# Make random colors repeatable: 54, 40,  4,99, 39,53,97 looks ok 
	print RCMD "barplot(x,beside=FALSE,col=sample(rainbow($NoOfcolours,start=$start, end=$end)),"; 
}
if (exists $parameters{"l"}){ print "with labels\n"; print RCMD "names.arg=data[,2],las=2,"; }	# print labels
print RCMD " border=NA,space=0,main='k=$k')\n";
#print RCMD " border=NA,space=0,main='k=$k seed=$seed')\n";
#}
close RCMD;

print "\n\nRunning R .. \n";
									# RUN R
`R CMD BATCH Rcmds.txt Rout.txt`;

print "\ndata:\t\t$datafile\nR output:\t$infile.Rout.txt\nR plots:\t$infile.pdf\n\n";
`mv Rcmds.txt $infile.Rcmds.txt`;
`mv Rout.txt $infile.Rout.txt`;


