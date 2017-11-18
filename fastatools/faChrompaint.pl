#!/usr/bin/perl


use warnings;
use strict;
use Getopt::Std;


my $program = 'faChrompaint.pl';			#name of script

my %parameters;						#input parameters
my ($infile,$cladefile,$colorfile,$ref, $infile_recog);
my $seq_recog = "a-z";				
my $W = 100000;						# window size
my $minW = $W*0.8;
my $exc = "'NCYC4146 1AA'";
my $NAcolor = "black";


my %ambcodes = (
	"M" => "AC",
	"R" => "AG",
	"W" => "AT",
	"S" => "CG",
	"Y" => "CT",
	"K" => "GT",
);

getopts('i:r:W:m:c:e:I:C:',\%parameters);

if (exists $parameters{"i"}) { $infile = $parameters{"i"}; }
if (exists $parameters{"I"}) { $infile_recog = $parameters{"I"}; }
if (exists $parameters{"r"}) { $ref = $parameters{"r"}; }
if (exists $parameters{"W"}) { $W = $parameters{"W"}; }
if (exists $parameters{"m"}) { $minW = $parameters{"m"}; }
if (exists $parameters{"c"}) { $cladefile = $parameters{"c"}; }
if (exists $parameters{"C"}) { $colorfile = $parameters{"C"}; }
if (exists $parameters{"e"}) { $exc = $parameters{"e"}; }

unless (((exists $parameters{"i"}) || (exists $parameters{"I"})) && (exists $parameters{"r"})) {
	print "\n USAGE: $program -i <infile> -r <ref>\n\n";
	print   "    -i\tfasta format infile\n";
	print   "    -I\tprefix for fasta format infile group (e.g. chr)\n";
	print   "    -r\treference seq name (e.g. P34048)\n";
	print   "    -W\tsliding window size [$W bp]\n";
	print   "    -m\tminimum good seq length for determining nearest strain [0.8*W]\n";
	print   "    -c\toptional cladefile (define clades to ID nearest clade)\n";
	print   "    -C\tfile of colors for clades (to go with -c)\n";
	print   "    -e\texclude strains matching a pattern [$exc]\n\n";
	print   " This script will summarise pairwise differences from a study strain in an alignment, and will optionally use R to paint chromosomes according to similarity with strains from known clades\n\n";
	print   " NOTE: ambiguity codes are treated as both bases (e.g. Y=C,Y=T)\n";
	print   " Format for cladefile: list of 'seqname cladename' (e.g. 'P87	4')\n\n";
	exit;
}

my $pdf = "$ref.paintchr.pdf";
my $outfile = "$ref.pDiffs.tsv";
open OUT, ">$outfile" or die "couldn't open $outfile : $!";

my (%clades,%colors,%multiples);	# need to handle multiple clade hits with smaller rectangles? %multiples: not written yet
if (defined $cladefile) {
	open CLADES, "<$cladefile" or die "couldn't open $cladefile : $!";
	while (<CLADES>) { 
		if (/^(\S+)\s+(\S+)/) { $clades{$1} = $2; }
	}
	open COLORS, "$colorfile" or die "\ncouldn't open file of colors ($colorfile) and need this for the -c option : $!\n\n";
	while (<COLORS>) { 
		if (/^(\S+)\s+(\S+)/) { $colors{$1} = $2; }
	}

}

my @infiles;
if (defined $infile) { push(@infiles,$infile); }
else { 
	opendir (DIR, ".") or die "couldn't open . : $!";

	my (%strains,@dirs);	
	while (defined (my $file = readdir(DIR))) {
		if ($file =~ /^$infile_recog/m) { push(@infiles,$file); }
	}
}
print "\nFiles to analyse: @infiles\n";

my $nearestout = "$ref.nearest.tsv";

open NEAREST, ">$nearestout" or die "couldn't open $nearestout : $!"; 
if (defined $cladefile) { print NEAREST "ref\trefclade\tinfile\twinpos\tnearestclade\tcladecolor\tneareststrains\n"; }
else { print NEAREST "ref\tinfile\twinpos\tneareststrains\n"; }


foreach my $infile (@infiles) {


	my ($seq_ref,@namesinfile,%seq,$length);
	($length,$seq_ref,@namesinfile) = fasta2strings($infile,$seq_recog);	# FASTA FORMAT SEQ
	%seq = %$seq_ref;

	my (@temp,%exc);
	if (defined $exc) {

		foreach my $name (@namesinfile) { 
			if ($exc =~ /\S+\s+.*?/) {				# if there is > 1 pattern match to make (INCOMPLETE)
				while ($exc =~ /(\S+)/g) { if ($name =~ /$1/) { $exc{$name}++; print "exclude $name\n"; } }
			}
			else { unless ($name =~ /$exc/) { push (@temp,$name); } }
		}
		unless (defined $temp[0]) { 
			foreach my $name (@namesinfile) { unless (defined $exc{$name}) { push (@temp,$name); } }
		}
		if (defined $temp[0]) { @namesinfile = @temp; print "$infile strains to study: @namesinfile\n"; }

	}

	unless (defined $seq{$ref}) { die "error: could not find reference sequence called $ref\n"; }

	if (defined $cladefile) { print OUT "file\tref\trefClade\tstrain\tstrainClade\tpos\tdiffs\tlengthNotN\tpDiff\n"; }
	else { print OUT "file\tref\tstrain\tpos\tdiffs\tlengthNotN\tpDiff\n"; }

	my (%minpdiff,%nearstrain);

	if ((defined $cladefile) && (!defined $clades{$ref})) { $clades{$ref} = "NA"; }

								 
	foreach my $name (@namesinfile) { 

		if ((defined $cladefile) && (!defined $clades{$name})) { $clades{$name} = "NA"; }

		if ($name eq $ref) { next; }					# skip the ref (compare other seqs to ref only)
		my $Wi = 0; my $diff = 0; my $T = 0; my $l = 0;				# position in window
		for (my $i=0; $i < @{$seq{$name}}; $i++) {
			if ((!defined $seq{$name}[$i]) || (!defined $seq{$ref}[$i])) { next; }	# avoid warnings with uneven chr lengths
			if (($seq{$name}[$i] eq "N") || ($seq{$name}[$i] eq "n") || ($seq{$ref}[$i] eq "N") || ($seq{$ref}[$i] eq "n")) { $Wi++; }	# skip low quality sequence in length
			elsif ($seq{$name}[$i] eq $seq{$ref}[$i]) { $Wi++; $l++; }	# no difference from reference
			elsif (($seq{$name}[$i] =~ /[MRWSYK]/) && ($seq{$ref}[$i] =~ /[$ambcodes{$seq{$name}[$i]}]/)) {	# seq is heterozygous and matches one allele of ref 
		#	print "seq: $seq{$name}[$i]\t$seq{$ref}[$i] = not a difference\n";
			}
			elsif (($seq{$ref}[$i] =~ /[MRWSYK]/) && ($seq{$name}[$i] =~ /[$ambcodes{$seq{$ref}[$i]}]/)) {	# seq is heterozygous and matches one allele of ref 
		#	print "seq: $seq{$name}[$i]\t$seq{$ref}[$i] = not a difference\n";
			}
			else { $diff++; $Wi++; $l++; }
			if ($Wi	== $W) {					# reached the end of a window	
				$T += $Wi;					# keep track of the true position
				if ($l == 0) { 
					if (defined $cladefile) { print OUT "$infile\t$ref\t$clades{$ref}\t$name\t$clades{$name}\t$T\t$diff\t$l\tNA\n"; }
					else { print OUT  "$infile\t$ref\t$name\t$T\t$diff\t$l\tNA\n"; }
					if (!defined $nearstrain{$T}) { $nearstrain{$T} = "NA"; } 
				}
				else {	
					my $pDiff = $diff/$l;
					if ($l >= $minW) {
						if ((defined $nearstrain{$T}) && ($nearstrain{$T} eq "NA")) { $minpdiff{$T} = $pDiff; $nearstrain{$T} = $name; }
						if (!defined $minpdiff{$T}) { $minpdiff{$T} = $pDiff; $nearstrain{$T} = $name; }
						elsif ($pDiff < $minpdiff{$T}) { $minpdiff{$T} = $pDiff; $nearstrain{$T} = $name; }
						elsif ($pDiff == $minpdiff{$T}) { $nearstrain{$T} .= ",$name" ; }
					} 
			#		else { $nearstrain{$T} = "NA"; }

					if (defined $cladefile) { print OUT "$infile\t$ref\t$clades{$ref}\t$name\t$clades{$name}\t$T\t$diff\t$l\t$pDiff\n"; }
					else { print OUT "$infile\t$ref\t$name\t$T\t$diff\t$l\t$pDiff\n"; }

				}
				$diff = 0; $Wi = 0; $l = 0;
			}
		}	
	}

	print "Summary of pairwise divergences between $ref and ".@namesinfile." sequences is in $outfile. The sliding window size was $W bp.\n\nThe sequence names compared to $ref were:\n@namesinfile\n\n";

# print out the most closely related strains for each sliding window position
# add the closest clade if a file of clades is included


	print "ref\trefclade\tinfile\twinpos\tnearestclade\tneareststrains\n";


	foreach my $T (sort { $a <=> $b } keys %nearstrain) {
	
		if (defined $cladefile) {					# USING USER-PROVIDED CLADE DEFINITIONS AND COLORS
			if ($nearstrain{$T} ne "NA") {
				my $nearclade;
				while ($nearstrain{$T} =~ /,?(\w+)/ig) { 
					if (defined $clades{$1}) {
						if (!defined $nearclade) { $nearclade = $clades{$1}; } 
						elsif ($nearclade ne $clades{$1}) { $nearclade = "multiple"; } 
					}
					else { warn "not in $cladefile: ($1)\n"; }
				}

				print "$ref\t$clades{$ref}\t$infile\t$T\t$nearclade\t$nearstrain{$T}\n"; 
				if (!defined $colors{$nearclade}) { 
					warn "No color specified for $nearstrain{$T} in Clade $nearclade, will use black\n";  
					$colors{$nearclade} = "black"; 
				}
				print NEAREST "$ref\t$clades{$ref}\t$infile\t$T\t$nearclade\t\"$colors{$nearclade}\"\t$nearstrain{$T}\n"; 
	
			} 							
			else { print "$ref\t$clades{$ref}\t$infile\t$T\tNA\t$nearstrain{$T}\n"; print NEAREST "$ref\t$clades{$ref}\t$infile\t$T\tNA\t$NAcolor\t$nearstrain{$T}\n"; }

		}
										# NO CLADE OR COLOR FILE
		else { 	print "$ref\t$infile\t$T\t$nearstrain{$T}\n"; print NEAREST "$ref\t$infile\t$T\t$nearstrain{$T}\n"; }

	}

	print "Summary of nearest strain(s) to $ref out of ".@namesinfile." sequences is in $nearestout. The sliding window size was $W bp.\n\nThe sequence names compared to $ref were:\n@namesinfile\n\n";

}

close OUT;
close NEAREST;



my $rcmdfile = "nearest$ref.R";
open RCMD, ">$rcmdfile" or die "couldn't open $rcmdfile : $!";
print RCMD "
rm(list=ls())
pdf(\"$pdf\")
data<-read.table(\"$nearestout\",header=T)
attach(data)
head(data)
tail(data)

par(mfrow=c(8,1),mar=c(1,1,1,1),yaxt='n',xlab='',bty=\"n\")

";


foreach my $infile (sort @infiles) {
	print RCMD "
plot(c(0,max(winpos)),c(0,100),type=\"n\",main=\"$ref.$infile\",yaxt=\'n\',xaxt=\'n\')
rect(winpos[infile==\"$infile\"]-100000,0,winpos[infile==\"$infile\"],100,col=as.vector(cladecolor[infile==\"$infile\"]))
";
}
close RCMD;

`R CMD BATCH --no-restore nearest$ref.R`;

					#################
					## SUBROUTINES ##
					#################	

				   #  IN ORDER OF APPEARANCE  #


## READ IN ALIGNMENT FROM FASTA FORMAT 
#
# NOTE SEQ IS RETURNED AS AN ARRAY, SLOW FOR LARGE DATA SETS
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
					#$seq{$name} = $seq;	# STRINGS
					
					@{$seq{$name}} = split (//, $seq); 	# ARRAYS
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

		#$seq{$name} = $seq;				# STRINGS

		@{$seq{$name}} = split (//, $seq); 		# ARRAYS		
	}
	else { print "WARNING: No sequence for $name in fasta file. It is not included in the analysis.\n"; }


	print "\nFound ".@names." sequences in fasta format in $datafile. ";
	if (defined $alignment) { print "Sequences are different lengths\n\n"; }
	else { print "All sequences are $length bp long.\n\n"; }

	return ($length,\%seq,@names);
}


