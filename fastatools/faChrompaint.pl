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
my $exc = "'NCYC4146 1AA SC5314_A'";
my $divcolor = "white";
my $NAcolor = "black";
my $warningcolor = "yellow";
my $maxintracladediffs = 1;		# if nearest distance is above this threshold, then nearest clade = "NA" and $NAcolor will be used
#my $maxintracladediffs = 0.001183812;			# if nearest distance is above this threshold, then nearest clade = "NA" and $NAcolor will be used


my %ambcodes = (
	"M" => "AC",
	"R" => "AG",
	"W" => "AT",
	"S" => "CG",
	"Y" => "CT",
	"K" => "GT",
);

getopts('i:r:W:m:c:e:I:C:M:w:p:',\%parameters);

if (exists $parameters{"i"}) { $infile = $parameters{"i"}; }
if (exists $parameters{"I"}) { $infile_recog = $parameters{"I"}; }
if (exists $parameters{"r"}) { $ref = $parameters{"r"}; }
if (exists $parameters{"W"}) { $W = $parameters{"W"}; $minW = $W*0.8; }
if (exists $parameters{"m"}) { $minW = $parameters{"m"}; }
if (exists $parameters{"c"}) { $cladefile = $parameters{"c"}; }
if (exists $parameters{"C"}) { $colorfile = $parameters{"C"}; }
if (exists $parameters{"e"}) { $exc = $parameters{"e"}; print "exclude: '$exc' if present\n"; }
if (exists $parameters{"M"}) { $maxintracladediffs = $parameters{"M"}; }
if (exists $parameters{"w"}) { $warningcolor = $parameters{"w"}; }


unless ((exists $parameters{"p"}) || (((exists $parameters{"i"}) || (exists $parameters{"I"})) && (exists $parameters{"r"}))) {
	print "\n USAGE: $program -i <infile> -r <ref>\n\n";
	print "    or: $program -I <prefix> -r <ref>\n\n";
	print "    or: $program -p <*.nearest.tsv> -c <cladefile> -C <colorfile>\n\n";


	print   "    -i\tfasta format infile\n";
	print   "    -I\tprefix for fasta format infile group (e.g. chr)\n";
	print   "    -p\t*.nearest.tsv file for Fast replotting with new colors\n";
	print   "    -r\tname of seq to compare to others (e.g. P34048)\n";
	print   "    -W\tsliding window size [$W bp]\n";
	print   "    -m\tminimum good seq length for determining nearest strain [0.8*W]\n";
	print   "    -c\toptional cladefile (define clades to ID nearest clade)\n";
	print   "    -C\tfile of colors for clades (to go with -c)\n";
	print   "    -M\tmaximum intraclade pDiffs [$maxintracladediffs]\n";
	print   "    -e\texclude strains matching a pattern [$exc]\n";
	print   "    -w\twarning color if nearest strain has no clade or color [$warningcolor]\n\n";
	print   " This script will:\n";
	print 	"    - summarise pairwise differences from a study strain in an alignment, \n";
	print	"    - optionally use R to paint chromosomes according to similarity with strains from known clades\n\n";
	print   " NOTE: ambiguity codes are treated as both bases (e.g. Y=C,Y=T)\n\n";
	print   " Format for cladefile: list of 'seqname cladename' (e.g. 'P87	4')\n";
	print   " Format for colorfile: list of 'cladename color' (e.g.'1 #dd1c77')\n\n";
	print   " Black = sequence quality is too low\n";
	print   " $warningcolor = nearest strain has no color in colorfile\n";
	print   " White = sequence is too diverged from all other sequences\n\n";
	print   " Without cladefile and colorfile you will only get a summary of closest sequences (no pictures)\n\n";
	exit;
}


my (%clades,%colors,%multiples,%multipleTs);	# need to handle multiple clade hits with smaller rectangles? %multiples: not written yet
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


my ($nearestout,$pdf,@infiles);

# USE A PAST .nearest.tsv FILE AS INPUT TO AVOID SLOW REPLOTTING WHEN ONLY CHANGING CLADE COLORS OR CLADE ASSIGNMENTS
# NOTE: you do need to re-run from scratch if you plan to change which strains you are excluding

if (defined $parameters{"p"}) { 
	$nearestout = $parameters{"p"};
	my %seenInfile;

	open IN, "<$nearestout" or die "couldn't open $nearestout : $!";

	my ($newnearestout,$header);
	
	while (<IN>) {
		if (/^(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)/m) {
			if ($1 eq "ref") { $header = $_; next; }		# header
			if ((defined $ref) && ($ref ne $1)) { die "error: there is more than 1 reference sequence in $nearestout\n"; }
			elsif (!defined $ref) { 				# first data line: identify reference sequence name
				$ref = $1;
				$newnearestout = "$ref.newnearest.tsv";
				open NEWNEAR, ">$newnearestout" or die "couldn't open $newnearestout : $!"; 
	 			print NEWNEAR $header; 
			}
			$ref = $1; my $infile = $2; my $T=$3; my $nstrains=$4; my $pDiff=$5;				# data

			unless ($nstrains =~ /,/) {									# only 1 near sequence
				if ((defined $clades{$nstrains}) && (defined $colors{$clades{$nstrains}})) {
					print NEWNEAR "$ref\t$clades{$ref}\t$infile\t$T\t$clades{$nstrains}\t\"$colors{$clades{$nstrains}}\"\t$nstrains\t$pDiff\n";
				}
				elsif ($nstrains eq "diverged") { 
					print NEWNEAR "$ref\t$clades{$ref}\t$infile\t$T\tdiverged\t\"$divcolor\"\t$nstrains\t$pDiff\n";
				}
				else { warn "unrecognized clade ($clades{$nstrains}) or strain ($nstrains)\n"; }
			}
			else {	
				my (@snames,%cladecount);
				while ($nstrains =~ /(\w+)/g) { push(@snames,$1); $cladecount{$clades{$1}}++; }
				my @cladenames = keys %cladecount;
				if (@cladenames == 1) { 									# similar strains are all from the same clade

					if (defined $colors{$cladenames[0]}) {

						print NEWNEAR "$ref\t$clades{$ref}\t$infile\t$T\t$cladenames[0]\t\"$colors{$cladenames[0]}\"\t$nstrains\t$pDiff\n";
					}
					else { warn "unrecognized clade ($cladenames[0]) or strain ($snames[0])\n"; }

				}
				else {												# there are multiple clades for this window
					push ( @{$multipleTs{"$infile"}}, $T );
					@{$multiples{"$infile$T"}} = @cladenames;

					print NEWNEAR "$ref\t$clades{$ref}\t$infile\t$T\tmultiple\t\"$colors{$cladenames[0]}\"\t$nstrains\t$pDiff\n";
				}
			}
			unless ($seenInfile{$infile}++) { push(@infiles,$infile); }	
		}
	}

	if ($nearestout =~ /(\S+?\.nearest)\.tsv/) { $pdf = "$1.pdf"; }
	else { $pdf = "$nearestout.pdf"; }

	close NEWNEAR;
	$nearestout = $newnearestout;
}

# OR ... figure out similarities from scratch

else {
if (defined $infile) { push(@infiles,$infile); }
elsif (defined $infile_recog) { 
	opendir (DIR, ".") or die "couldn't open . : $!";

	my (%strains,@dirs);	
	while (defined (my $file = readdir(DIR))) {
		if ($file =~ /^$infile_recog/m) { push(@infiles,$file); }
	}
}
print "\nFiles to analyse: @infiles\n";

my $outfile = "$ref.pDiffs.tsv";
open OUT, ">$outfile" or die "couldn't open $outfile : $!";

$nearestout = "$ref.nearest.tsv";
my $seenFilter;

open NEAREST, ">$nearestout" or die "couldn't open $nearestout : $!"; 
if (defined $cladefile) { print NEAREST "ref\trefclade\tinfile\twinpos\tnearestclade\tcladecolor\tneareststrains\tnearestpdiff\n"; }
else { print NEAREST "ref\tinfile\twinpos\tneareststrains\tnearestpdiff\n"; }

my $excpdf; my %seene;
foreach my $infile (sort @infiles) {


	my ($seq_ref,@namesinfile,%seq,$length);
	($length,$seq_ref,@namesinfile) = fasta2strings($infile,$seq_recog);	# FASTA FORMAT SEQ
	%seq = %$seq_ref;

	my (@temp,%exc);							# SEQUENCES ARE EXCLUDED WHEN INITIALLY READING THE ALIGNMENT (therefore cannot change strain strain exclusions after creating nearest.tsv)
	if (defined $exc) {

		foreach my $name (@namesinfile) { 
		#	print "do I need to exclude '$name'?\n";
			if ($exc =~ /\S+\s+.*?/) {				# if there is > 1 pattern match to make 
				while ($exc =~ /(\S+)/g) { 
				#	print "checking for '$1' to exclude\n";
 					my $e = $1;
					if ($name =~ /$e/) { 
						$exc{$name}++; print "exclude $name\n"; unless ($seene{$e}++) { $excpdf .= $e; } 
					} 
				}
			}
			else { unless ($name =~ /$exc/) { push (@temp,$name); } }
		}
		unless (defined $temp[0]) { 
			foreach my $name (@namesinfile) { unless (defined $exc{$name}) { push (@temp,$name); } }
		}
		if (defined $temp[0]) { @namesinfile = @temp; print "$infile strains to probe with: @namesinfile\n"; }

	}

	unless (defined $seq{$ref}) { die "error: could not find reference sequence called $ref\n"; }

	if (defined $cladefile) { print OUT "file\tref\trefClade\tstrain\tstrainClade\tpos\tdiffs\tlengthNotN\tpDiff\n"; }
	else { print OUT "file\tref\tstrain\tpos\tdiffs\tlengthNotN\tpDiff\n"; }

	my (%minpdiff,%nearstrain,%nearestpdiff,%seen);

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
					if (!defined $nearstrain{$T}) { $nearstrain{$T} = "NA"; $nearestpdiff{$T} = "NA"; } 
				}
				else {	
					my $pDiff = $diff/$l;
					if (($l >= $minW) && ($pDiff < $maxintracladediffs)) {
						if ((defined $nearstrain{$T}) && ($nearstrain{$T} eq "NA")) {
							unless ($seenFilter++) { print "filtering regions that have diverged over $maxintracladediffs from all other sequences\n"; } 
							$minpdiff{$T} = $pDiff; $nearstrain{$T} = $name; $nearestpdiff{$T} = $pDiff;
						}
						if (!defined $minpdiff{$T}) { $minpdiff{$T} = $pDiff; $nearstrain{$T} = $name; $nearestpdiff{$T} = $pDiff; }
						elsif ($pDiff < $minpdiff{$T}) { $minpdiff{$T} = $pDiff; $nearstrain{$T} = $name; $nearestpdiff{$T} = $pDiff; }
						elsif ($pDiff == $minpdiff{$T}) { $nearstrain{$T} .= ",$name" ; }
					} 
					elsif (($pDiff > $maxintracladediffs) &&  (!defined $nearstrain{$T})) { 
						$nearstrain{$T} = "diverged"; $nearestpdiff{$T} = $pDiff; 
					}
					elsif (!defined $nearstrain{$T}) { $nearstrain{$T} = "NA"; $nearestpdiff{$T} = $pDiff; }

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


	print "ref\trefclade\tinfile\twinpos\tnearestclade\tneareststrains\tnearestpdiff\n";


	foreach my $T (sort { $a <=> $b } keys %nearstrain) {
	
		if (defined $cladefile) {					# USING USER-PROVIDED CLADE DEFINITIONS AND COLORS
			if ($nearstrain{$T} eq "diverged") {
				print "$ref\t$clades{$ref}\t$infile\t$T\tdiverged\t$nearstrain{$T}\t$nearestpdiff{$T}\n"; print NEAREST "$ref\t$clades{$ref}\t$infile\t$T\tdiverged\t$divcolor\t$nearstrain{$T}\t$nearestpdiff{$T}\n";
			}
			elsif ($nearstrain{$T} ne "NA") {
				my $nearclade;
				while ($nearstrain{$T} =~ /,?([0-9a-z\.]+)/ig) { 
					if (defined $clades{$1}) {
						if (!defined $nearclade) { $nearclade = $clades{$1}; } 
						elsif ($nearclade ne $clades{$1}) { 
							if (!defined $multiples{"$infile$T"}) { 
								push ( @{$multipleTs{"$infile"}}, $T ); 
								push ( @{$multiples{"$infile$T"}}, $nearclade ); $seen{"$infile$T$nearclade"}++;
								push ( @{$multiples{"$infile$T"}}, $clades{$1} ); $seen{"$infile$T$clades{$1}"}++;
 
								$nearclade = "multiple";
							}							
							else { 
								unless ($seen{"$infile$T$clades{$1}"}++) {
									push (@{$multiples{"$infile$T"}},  $clades{$1}); 
								}
							}
						} 
					}
					else { warn "not in $cladefile: ($1)\n"; }
				}

				print "$ref\t$clades{$ref}\t$infile\t$T\t$nearclade\t$nearstrain{$T}\t$nearestpdiff{$T}\n"; 
				if (!defined $colors{$nearclade}) { 
					warn "No color specified for $nearstrain{$T} in Clade $nearclade, will use $warningcolor\n";  # NO-COLOR-SPECIFIFED WARNING SHOULD BE DIFFERENT THAN LOW QUALITY COLOR 
					$colors{$nearclade} = $warningcolor; 
				}
				print NEAREST "$ref\t$clades{$ref}\t$infile\t$T\t$nearclade\t\"$colors{$nearclade}\"\t$nearstrain{$T}\t$nearestpdiff{$T}\n"; 
	
			} 							
			else { print "$ref\t$clades{$ref}\t$infile\t$T\tNA\t$nearstrain{$T}\t$nearestpdiff{$T}\n"; print NEAREST "$ref\t$clades{$ref}\t$infile\t$T\tNA\t$NAcolor\t$nearstrain{$T}\t$nearestpdiff{$T}\n"; }

		}
										# NO CLADE OR COLOR FILE
		else { 	print "$ref\t$infile\t$T\t$nearstrain{$T}\t$nearestpdiff{$T}\n"; print NEAREST "$ref\t$infile\t$T\t$nearstrain{$T}\t$nearestpdiff{$T}\n"; }

	}

	print "Summary of nearest strain(s) to $ref out of ".@namesinfile." sequences is in $nearestout. The sliding window size was $W bp.\n\nThe sequence names compared to $ref were:\n@namesinfile\n\n";

}

close OUT;
close NEAREST;


# GIVE THE OUTPUT PDF A SENSIBLE NAME THAT KEEPS TRACK OF THE OPTIONS SELECTED
$pdf = $ref;
if (defined $parameters{"M"}) { $pdf .= "paintchrM"; }
else { $pdf .= "paintchr"; }
if ((defined $exc) && (defined $excpdf)) { $pdf .= "e$excpdf"; }
$pdf .= ".pdf";

}






# PRINT R CMDS FOR VISUALIZING THE NEAREST SEQ IN IN EVERY WINDOW IN THE GENOME
my $rcmdfile = "nearest$ref.R";
open RCMD, ">$rcmdfile" or die "couldn't open $rcmdfile : $!";
print RCMD "
rm(list=ls())
pdf(\"$pdf\")
data<-read.table(\"$nearestout\",header=T)
attach(data)
head(data)
tail(data)

par(mfrow=c(8,1),mar=c(1,1,1,1),yaxt='n',bty=\"n\")

";


foreach my $infile (sort @infiles) {

# PRINT STRAIGHTFORWARD GENOME PLOTS USING USER-SPECIFIED COLOUR FOR BINS WITH TIED MATCHES

	print RCMD "

# PRINT STRAIGHTFORWARD GENOME PLOTS USING USER-SPECIFIED COLOUR FOR BINS WITH TIED MATCHES
# data from $nearestout
plot(c(0,max(winpos)),c(0,100),type=\"n\",yaxt=\'n\',xaxt=\'n\',main=\"$ref.$infile\")
rect(winpos[infile==\"$infile\"]-$W,0,winpos[infile==\"$infile\"],100,col=as.vector(cladecolor[infile==\"$infile\"]))
";


# OVERWRITE WITH SMALLER RECTANGLES TO SHOW MULTIPLE COLORS IN CASE OF TIES 

	if (defined $multipleTs{"$infile"}[0]) { 
		my @windows = @{$multipleTs{"$infile"}};
		print "OK $infile contains equally similar sequences at ".@windows." windows: @windows\n"; 
		foreach my $T (@windows) {
			my @clades = sort @{$multiples{"$infile$T"}};
			print "$infile $T [@clades]\n";
			my $n = @clades;
			my $low = 0;
			my $high = 100/$n;
			for (my $i=0; $i < @clades; $i++) { 
				$low += 0+($i*100/$n);
				$high += $i*100/$n;
				print RCMD "

# OVERWRITE WITH SMALLER RECTANGLES TO SHOW MULTIPLE COLORS IN CASE OF TIES 
# data is added in using perl and is not contained in $nearestout
rect($T-$W,$low,$T,$high,col=\"$colors{$clades[$i]}\")
";
			}
		}
	}
}
close RCMD;


# RUN R

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


