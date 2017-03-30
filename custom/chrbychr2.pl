#!/usr/bin/perl


use warnings;
use strict;

my $infile = "Rout.txt";
my $bsmin = 7000;			# only consider clusters with more than this number of bootstraps
my $outfile = "chrbychr$bsmin.txt";
my $Rcmds = "Rcmds.txt";
my $Rout = "Rout2.txt";
my $habitatFile = "habitatSimple.txt";


				# the popns defined in Liti et al 2009 + subpopulations for European and USA
				# write code to handle subpopulations properly: 
				# 	- European or USA when subpopulations not defined
my %popns = (
	"YJM978m" => "WINE",
	"NCYC110m" => "WA",
	"YPS606m" => "PAO",
	"Y9m" => "sake",
	"UWOPS05_2173m" => "MAL",
	"SDO1s1" => "NCO",
	"ZP560s1" => "EUO",
);
my %subpopns = (
	"WINE" => "European",
	"EUO" => "European",
	"NCO" => "USA",
	"PAO" => "USA",
);

my %habitat;					# read in summary of habitats for Fisher's Exact Test
open HABITAT, "$habitatFile" or die "couldn't open $habitatFile : $!";
while (<HABITAT>) { 
	if (/^(\S+)\s+(\S+)/) {
		my $strain = $1;
		my $habitat = $2;
		if ($habitat ne "HabitatSimple") { $habitat{$strain} = $habitat; }
	}
	else { print "unrecognised line in $habitatFile : $_"; } 
}



open RESULTS, "<$infile" or die "couldn't open $infile : $!";

my ($locus, %strain, $prevCluster, $prevBs, %longestPopn,%longSubPop, %allpopns, @loci);
while (<RESULTS>) {
	if (/prop\.part\(cen(\d+)bs\$trees\)/) { 
		if (%strain) {						# if this is not the first locus
									# summarise results from previous locus
			my $allpopnsref = strainSummary($locus,\%strain,\%longestPopn,\%longSubPop,\%allpopns);	
			%allpopns = %$allpopnsref;
			%longestPopn = ();
			%longSubPop = ();
			%strain = ();
		}

		$locus = $1;						# starting a new chromosome analysis
		push(@loci, $locus);	
	}
	if ((defined $locus) && (/(\d+):\s+cen$locus\_(\S+)/)) { 
		$strain{$1} = $2;
	}
	if (/^==> (\d+) time\(s\)\s*:\s*\[\d+\]\s*(.*)/m) {		# reading 1st line of bootstraps
		if (defined $prevBs) {		
			printCluster($prevBs,$prevCluster);
		}
		$prevBs = $1;
		$prevCluster = $2;
	}
	elsif (/^\[\d+\]\s*(.*)/m) { $prevCluster .= " $1"; }		# reading subsequent line of bootstraps

}
if (defined $prevBs) {		my $longestPopnRef = 
	printCluster($prevBs,$prevCluster);
}
my $allpopnsref = strainSummary($locus, \%strain,\%longestPopn,\%longSubPop,\%allpopns);	# summarise results from final locus 
%allpopns = %$allpopnsref;

open OUT, ">$outfile" or die "couldn't open $outfile : $!";
print OUT "Strain\tHabitat\tUsedAsRef\t";
foreach my $chr (@loci) { print OUT "chr$chr\t"; };
print OUT "chrbychr70\tchrbychr70details\n";

foreach my $strain (sort values %strain) {
	print OUT "$strain\t$habitat{$strain}\t";
	if (defined $popns{$strain}) { print OUT "ref\t"; }
	else { print OUT "nonref\t"; }
	
#	print "$strain\t@{allpopns{$strain}}";
	my (%seen, @uniqpopns);
	foreach my $chr (@loci) {
		if (!defined $allpopns{$strain}[$chr]) { print OUT "na\t"; }	# "na" if missing data
		else {
		       	unless ($seen{$allpopns{$strain}[$chr]}++) { push(@uniqpopns, $allpopns{$strain}[$chr]); }	
			print OUT "$allpopns{$strain}[$chr]\t"; 
		}			# print results
	}
	if (@uniqpopns > 1) {						# more than 1 designation (but that includes undef or European/USA)
		my $P;
		foreach my $p (@uniqpopns) {
			if (!defined $P) { $P = $p; }				# the first popn descriptor
			elsif (($P ne $p) && ($P eq 'undef')) { $P = $p; }	# replace 'undef' with a definition if available
			elsif (($P ne $p) && ($p ne 'undef')) {			# if > 1 popn name is defined 
				if ((exists $subpopns{$p}) && ($subpopns{$p} eq $P)) { 
					$P = $p; 			# use subpopn name if General and Subpopn name are available
				}		
				elsif ((exists $subpopns{$P}) && ($subpopns{$P} eq $p)) { 
					#no change			# use subpopn name if General and Subpopn name are available
				} 
				else { $P = "admixed"; }
			}
		}
		my $summary = join ("/",@uniqpopns);
		print OUT "$P\t$summary"; 
	}
	elsif (@uniqpopns == 1) { print OUT "$uniqpopns[0]\t$uniqpopns[0]"; }
       	else { print "error: $strain <@uniqpopns> [@{$allpopns{$strain}}]\n"; } 	
	print OUT "\n";
}
print "\n\nfinal summary of chr-by-chr population assignments for every strain is in $outfile\n\n";

#exit;

# MAKE TABLES AND RUN FISHER'S EXACT TESTS IN R

open RCMD, ">$Rcmds" or die "couldn't open $Rcmds : $!";
print RCMD "
data<-read.table(\"$outfile\",header=T)
attach(data)
head(data)


# Fisher's exact test on populations (USA, European etc) vs admixed/non-admixed
h<-chrbychr70
levels(h)
levels(h)<-c('admixed', rep('nonadmixed',7))
t<-table(Habitat,h)
t
fisher.test(t)


# Excluding oaks: Fisher's exact test on populations (USA, European etc) vs admixed/non-admixed
t2<-t[c(1:3,5),]
t2
fisher.test(t2)

";
close RCMD;

print "\n\nRunning R .. \n";
									# RUN R
`R CMD BATCH $Rcmds $Rout`;

print "\ndata:\t\t$outfile \nR output:\t$Rout\n\n";

system("cat $Rout"); 



					#################
					## SUBROUTINES ##
					#################	

				   #  IN ORDER OF APPEARANCE  #


## READ CLUSTER, REPORT THOSE WITH HIGH BOOTSTRAPS,  convert clusters of numbers into clusters of strains
#  report both sides of the bipartition, show clean popns from Liti et al 2009

sub printCluster {
	       	my $bs = shift;
		my $clusters = shift;	
		if ($bs >= $bsmin) { 	# if bootstrap support is above $bsmin	
					# convert clusters of numbers into clusters of strains
			my @clusters = split (/\s+/, $clusters);
		#      	print "\nA\t$locus\t$bs\t";	 # comment this out when not troubleshooting
	

# report both sides of the bipartition
# The A side
			my %seen;	
			my %seenPopnA;	# is this a diagnostic cluster?
			my @strainsA;
			foreach my $no (@clusters) {
				unless (defined $strain{$no}) { next; }
			       	my $strain = $strain{$no};	# convert numbers to strains
				push(@strainsA, $strain);	# list of strains in this cluster
				$seen{$strain}++;		# track which strains are in this side of the bipartition
		#		print "$strain\t"; 	# comment this out when not troubleshooting

				if (exists $popns{$strain}) {		# reference strain 
					$seenPopnA{$popns{$strain}}++;
				}
				else { 
				#	print "$strain "; 
				}
			}
			my @popns = keys %seenPopnA;		# cluster is from only 1 popn?
			my $defined = @popns;		# cluster is from only 1 popn?
		#	print "\nNo. of popns=$defined [@popns]\n"; 	# comment this out when not troubleshooting

		
		
			if ($defined == 1) { 
				my $longRef = diagnostic($popns[0],\%longestPopn,@strainsA);
				%longestPopn = %$longRef;
			}					# or 2 subpopulations from same popn
			elsif  (($defined == 2) && (exists $subpopns{$popns[0]}) && (exists $subpopns{$popns[1]}) && ($subpopns{$popns[0]} eq $subpopns{$popns[1]})) {

				my $longRef = diagnostic($subpopns{$popns[0]},\%longSubPop,@strainsA);
				%longSubPop = %$longRef;
			}

# report both sides of the bipartition
# The B side
			my %seenPopnB;				# is this a diagnostic cluster?
			my @strainsB;
		#	print "\nB\t$locus\t$bs\t"; # comment this out when not troubleshooting
			foreach my $strain (values %strain) {	
				unless (defined $seen{$strain}) {
				       	push(@strainsB, $strain);	# list of strains in this cluster
		#			print "$strain\t"; 	# comment this out when not troubleshooting
	
					if (exists $popns{$strain}) { 
						$seenPopnB{$popns{$strain}}++;	
					}
					else { 
					}
				}
			}
			@popns = keys %seenPopnB;		# cluster is from only 1 popn?
			if (@popns == 1) { 
				my $longRef = diagnostic($popns[0],\%longestPopn,@strainsB);
				%longestPopn = %$longRef;

			}					# or 2 subpopulations from same popn
			elsif  ((@popns == 2) && (exists $subpopns{$popns[0]}) && (exists $subpopns{$popns[1]}) && ($subpopns{$popns[0]} eq $subpopns{$popns[1]})) {

				my $longRef = diagnostic($subpopns{$popns[0]},\%longSubPop,@strainsB);
				%longSubPop = %$longRef;
			}
	#	print "END cluster\n";	# comment this out when not troubleshooting
		} 

}

sub diagnostic {

		my $pop = shift;
		my $longestPopnRef = shift;
		my @strains = @_;
		
		my %longestPopn = %$longestPopnRef;
		
		if (!defined $pop) { print "Pop is undefined: $pop \t @strains\n"; }
	#	print " DIAGNOSTIC ($pop)";
							# is this the longest diagnostic cluster?
		unless (defined @{$longestPopn{$pop}}) { @{$longestPopn{$pop}} = @strains; }
							# warn if > 1 cluster of same length for a popn (first cluster is chosen)
		elsif (@strains == @{$longestPopn{$pop}}) { 
			warn "\n\n> 1 cluster has same length for $pop. The first cluster is chosen.\n\@strains: @strains \@{\$longestPopn{$pop}}) : @{$longestPopn{$pop}})\n";  		} 
		elsif (@strains >  @{$longestPopn{$pop}}) {  @{$longestPopn{$pop}} = @strains; }


		return \%longestPopn;
}

sub strainSummary {
	my $chr = shift;
	my $strainRef = shift;
	my $longestPop = shift;
	my $longSubPop = shift;
	my $allpopnsref = shift;

	my @strains = sort values %$strainRef;
	my %longestPop = %$longestPop;
	my %longSubPop = %$longSubPop;
	my %allpopns = %$allpopnsref;

	my %seen;
	my $checkstrainNo = 0;

	print "\n\nSummary of strain diagnosis for $locus:\n";
	foreach my $pop (sort keys %longestPop) {
		foreach my $strain (@{$longestPop{$pop}}) { $seen{$strain}++; }	# have all strains been categorised?
	}

	foreach my $strain (@strains) {
		unless ($seen{$strain}) { 							# if uncategorised:
			foreach my $pop (sort keys %longSubPop) {				# define as European or USA
			#	print "hello: $pop [@{$longSubPop{$pop}}]\n";
				foreach my $s (@{$longSubPop{$pop}}) {
					if ((defined $s) && ($s eq $strain)) { 
						push (@{$longestPop{$pop}},$strain); 
						$seen{$strain}++; 
					}
				}
			}									# or
		       	unless ($seen{$strain}) {						# use Liti et al category
				if ($popns{$strain}) { push (@{$longestPop{$popns{$strain}}}, $strain); }
				else { push(@{$longestPop{"undef"}}, $strain); }		#  or "undef" 
			}
		}	
	}

	foreach my $pop (sort keys %longestPop) {
		$checkstrainNo += @{$longestPop{$pop}};
		print "$pop\t".@{$longestPop{$pop}}."\t";
		foreach my $strain (@{$longestPop{$pop}}) { 	
			print "$strain "; 
									# save a summary hash of arrays with every strain assignment
			$allpopns{$strain}[$chr] = $pop;
		}
		print "\n";
	}
	print "Total: $checkstrainNo\n\n";

	return \%allpopns;
}
