#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

my $program = 'vcf2allelePlot.pl';					#name of script

my %parameters;							#input parameters
my ($infile,$outfile,$gffile);
my $qual = 40;
my $RcmdFile = "vcf2allelePlot.Rcmds";
my $mwsize = 5000;

getopts('i:o:q:g:mh',\%parameters);

if (exists $parameters{"i"}) { $infile = $parameters{"i"}; }
if (exists $parameters{"q"}) { $qual = $parameters{"q"}; }
if (exists $parameters{"g"}) { $gffile = $parameters{"g"}; }
if (exists $parameters{"o"}) { $outfile = $parameters{"o"}; }
elsif (defined $infile) { 
	if ($infile =~ /^(.+)\.\S+?$/m) { $outfile = "$1.calls"; }
	else { $outfile = "$infile.calls"; }
} 

unless (exists $parameters{"i"}) {
	print "\n USAGE: $program -i '<vcf file>'\n\n";
	print   "    -i\tvcf file from mpileup -u + bcftools call -c\n";
	print   "    -q\tminimum phred-scaled quality [$qual]\n";
	print   "    -g\tgff file of annotations [none]\n";
	print   "    -m\tshow mode in 5kb non-overlapping sliding windows\n";
	print   "    -h\tshow heterozygosity in 5kb non-overlapping sliding windows\n";
	print 	"    -o\tprefix for outfiles [prefix of vcf file]\n\n";
	exit;
}

# READ THE GFF OF ANNOTATIONS IF THERE IS ONE

my (%atype,%astart,%aend,%seentype);		# chr will be the keys in these hashes of arrays
my (@chr, @types);
if (defined $gffile) { 
	open GFF, "<$gffile" or die "couldnt open $gffile : $!"; 
	print "\nReading in annotations from $gffile\n";
	while (<GFF>) {
		if (/^(\S+)\s+\S+\s+(\S+)\s+(\d+)\s+(\d+)/m) {
			$seentype{$2}++;
			push (@{$atype{$1}},$2);
			push (@{$astart{$1}},$3);
			push (@{$aend{$1}},$4);			
		}
		else { warn "\n**Unrecognized format for $gffile.**\n\n"; }
			
	}
	close GFF;
	@chr = sort keys %atype;
	@types = sort keys %seentype;
	print "   Found annotations for ".@chr." chromosomes: @chr\n";
	print "   There are ".@types." types of annotation: @types\n";

}


# Extract relevant information for the allele plot 
# and Read in information relevant for the R commands (e.g. which chromosomes are present)
#
open DATA, "<$infile" or die "couldn't open $infile : $!";

open OUT, ">$outfile" or die "couldn't open $outfile : $!";

print OUT "chr\tpos\tREF\tALT\tQUAL\tREFfwd\tREFrev\tALTfwd\tALTrev\tpALT\ttype\n";

print "\nReading and printing data from $infile to $outfile ..\n";

my %chr;
	
while (<DATA>) { 
								# find the lines with data 
	if (/^(\S+)\s+(\d+)\s+\.\s+(\S+)\s+(\S+)\s+(\S+)\s+\.\s+\S+?DP4=(\d+),(\d+),(\d+),(\d+)/m) {
		if ($4 eq ".") { next; }			# skip invariant sites
		my $pAlt = ($8+$9)/($6+$7+$8+$9);
		print OUT "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\t$pAlt\t";
		$chr{$1}++;					# a hash storing chromosome names and the number of variant sites for each
		my $ref = $3;
		my $alt = $4;
		my $variant = "snp";				# a non-SNP variant
		unless (($ref =~/^\S$/m) && ($ref =~ /^\S$/m)) { $variant = "indel"; }
		print OUT "$variant\n";
	}
}
close OUT;


# Print Rcmds to run with R CMD BATCH
open RCMD, ">$RcmdFile" or die "couldn't open $RcmdFile : $!";

print RCMD "
rm(list=ls())
data<-read.table(\"$outfile\",header=T)
attach(data)
head(data)

				# A function for estimating the mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

";
								# add an extra plot for heterozygosity if requested
if (defined $parameters{"h"}) { print RCMD "par(mfrow=c(3,1),cex=0.5)\npar(xpd=TRUE)\n";  } 
else { print RCMD "par(mfrow=c(2,1),cex=0.5)\n"; }

foreach my $chr (sort keys %chr) {

	print RCMD "

				# prepare a sliding window vector for estimating mode
x<-ceiling(max(pos[chr==\"$chr\"&QUAL>=40])/100000)
W<-1:((x/($mwsize/1000))*100)*$mwsize-($mwsize/2)
head(W)
tail(W)
modepALTsnp<-0
modepALTindel<-0
diff_snp<-0
err_snp<-0
het_snp<-0

for(i in 1:length(W)) { 
	modepALTsnp[i] <- Mode(pALT[pos>(W[i]-($mwsize/2))&pos<=(W[i]+($mwsize/2))&chr==\"$chr\"&QUAL>=$qual&type==\"snp\"]) 
	diff_snp[i] <- sum(pALT[pos>(W[i]-($mwsize/2))&pos<=(W[i]+($mwsize/2))&chr==\"$chr\"&QUAL>=$qual&type==\"snp\"]>0.8) 
	err_snp[i] <- sum(pALT[pos>(W[i]-($mwsize/2))&pos<=(W[i]+($mwsize/2))&chr==\"$chr\"&QUAL>=$qual&type==\"snp\"]<0.2) 
	het_snp[i] <- sum(pALT[pos>(W[i]-($mwsize/2))&pos<=(W[i]+($mwsize/2))&chr==\"$chr\"&QUAL>=$qual&type==\"snp\"])-err_snp[i]-diff_snp[i]
}
summary(modepALTsnp)
head(cbind(W,modepALTsnp))
tail(cbind(W,modepALTsnp))
head(cbind(W,diff_snp/$mwsize))
head(cbind(W,err_snp/$mwsize))

for(i in 1:length(W)) { modepALTindel[i] <- Mode(pALT[pos>(W[i]-($mwsize/2))&pos<=(W[i]+($mwsize/2))&chr==\"$chr\"&QUAL>=$qual&type==\"indel\"]) }


				# MAKE SEPARATE STACKED PLOTS FOR SNPs AND INDELS
				# SNPs
plot(c(0,max(W)),c(0,1),main=\"$infile: $chr freq of alternate alleles\",sub=\"SNPs Q$qual+\",xlab=\"position\",xaxt=\"n\",ylim=c(0,1),xlim=c(1,max(pos[chr==\"$chr\"])),type = \"n\")
points(pos[chr==\"$chr\"&QUAL>=$qual&type==\"snp\"],pALT[chr==\"$chr\"&QUAL>=40&type==\"snp\"],pch=20,col=\"black\")
axis(1, xaxp=c(0, signif(max(pos[chr==\"$chr\"]),3), 20))
abline(h=0.5)
abline(h=mean(pALT[chr==\"$chr\"&QUAL>=$qual&type==\"snp\"]),col=\"orange\")
\n";

	if (defined $parameters{'m'}) { showmode($chr,"snp"); }	# show mode if requested
	if (defined $parameters{'g'}) { annotate($chr); }	# show annotations on SNP plot if requested 
	
	if (defined $parameters{'h'}) { 			# make a new plot of H if requested
								# legend format for 3 plots per page
		print RCMD "legend(0,-0.16,c(round(mean(pALT[chr==\"$chr\"&QUAL>=$qual&type==\"snp\"]),3),0.5),lty=c(1,1),col=c(\"orange\",\"black\"),title=\"mean\",bty=\"n\")\n";
		showH($chr,"snp");
		annotate($chr); 				# show annotations on H plot if requested 
	}	
	else { print RCMD "legend(0,0.15,c(round(mean(pALT[chr==\"$chr\"&QUAL>=$qual&type==\"snp\"]),3),0.5),lty=c(1,1),col=c(\"orange\",\"black\"),title=\"mean\",bty=\"n\")\n"; }				# legend format for 2 plots per page		
	
	print RCMD "
				# INDELs
plot(c(0,max(W)),c(0,1),main=\"$infile: $chr freq of alternate alleles\",sub=\"Indels Q$qual+\",xlab=\"position\",xaxt=\"n\",ylim=c(0,1),xlim=c(1,max(pos[chr==\"$chr\"])),type = \"n\")	
points(pos[chr==\"$chr\"&QUAL>=$qual&type==\"indel\"],pALT[chr==\"$chr\"&QUAL>=40&type==\"indel\"],pch=20,col=\"black\")	
axis(1, xaxp=c(0, signif(max(pos[chr==\"$chr\"]),3), 20))
abline(h=0.5)
abline(h=mean(pALT[chr==\"$chr\"&QUAL>=$qual&type==\"indel\"]),col=\"orange\")
legend(0,-0.16,c(round(mean(pALT[chr==\"$chr\"&QUAL>=$qual&type==\"indel\"]),3),0.5),lty=c(1,1),col=c(\"orange\",\"black\"),title=\"mean\",bty=\"n\")
";

	if (defined $parameters{'m'}) { showmode($chr,"indel"); }	# show mode if requested
	if (defined $parameters{'g'}) { annotate($chr); }		# show annotations on INDEL plot if requested (in progress)
	
	print RCMD "rm(W,modepALTsnp,modepALTindel)\n";			# CLEAN UP


}

close RCMD;

# RUN THE R COMMANDS WITH OUTPUTS GOING OUT TO THE DEFAULT FILE NAMES

print "Running $RcmdFile commands in R ..\n";
`R CMD BATCH $RcmdFile`;
print "Done. R output is in $RcmdFile".".Rout and plots are in Rplots.pdf\n\n";


# SUBROUTINES 

sub annotate {
						# annotate the plot with slightly transparent colored rectangles
	my $chr = shift;

	print RCMD "par(xpd=F)\n";	# don't print annotations outside the plot	
	for (my $i=0; $i<@{$astart{$chr}}; $i++) {
		my $j;				# use a number from 1 to n for type color	
		for ($j=0; $j<@types; $j++) { if ($atype{$chr}[$i] eq $types[$j]) { $j += 1; last; } }	
			
		print RCMD "rect($astart{$chr}[$i],0,$aend{$chr}[$i],1,col=rainbow(".@types.",alpha=0.3)[$j],border=rainbow(".@types.",alpha=0.3)[$j])\n"; 	

	}
	print RCMD "par(xpd=T)\n";	# do print legend outside the plot	


	for (my $j=1; $j<=@types; $j++) {	
		print RCMD "legend(W[60],-0.16,legend=\"$types[$j-1]\",col=rainbow(".@types.",alpha=0.3)[$j],bty=\"n\")\n";
	}

}	

sub showmode {
	my $chr = shift;
	my $type = shift;

	print RCMD "lines(W,modepALT$type,col=\"green\")\n";
	print RCMD "legend(W[60],-0.16,Mode(pALT[chr==\"$chr\"&QUAL>=$qual&type==\"$type\"]),lty=1,col=\"green\",title=\"mode\",bty=\"n\")\n";

}

sub showH {
	my $chr = shift;
	my $type = shift;

	print RCMD "
my<-max((het_$type+diff_$type+err_$type)/$mwsize)
plot(c(0,max(W)),c(0,max((het_$type+diff_$type+err_$type)/$mwsize)),type='n',xlab=\"position\",main=\"$infile: $chr sliding window ($mwsize) of heterozygosity, homozygous diffs and error\",sub=\"$type Q$qual+\",xaxt=\"n\")
axis(1, xaxp=c(0, signif(max(pos[chr==\"$chr\"]),3), 20))
\n";
	print RCMD "lines(W,diff_$type/$mwsize,col=\"blue\")\n";
	print RCMD "lines(W,het_$type/$mwsize,col=\"red\")\n";
	print RCMD "lines(W,err_$type/$mwsize,col=\"grey\")\n";

	print RCMD "legend(W[60],my,c(\"diff_$type\",\"het_$type\",\"err_$type\"),lty=c(1,1,1),col=c(\"blue\",\"red\",\"grey\"),bty=\"n\")\n";

}				

