#!/usr/bin/perl

use warnings;
use strict;

my $program = 'runFinal2.pl';					#name of script

my @chr = qw(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chrR);
my $outname = "roparsSuppTable3";
my $foutname = "roparsFiltered.mfa";

my $dir = "/lustre1/doudab/candidAnalysis/final2";
my $bootstraps = 100;
my $ncbireaddesc = "ncbisra_result.tsv";
my $ebireaddesc = "../Ropars_etal18/PRJNA432884.txt";
my $clades = "clades";
#my $withinCladeMax = 0.0006604028;
#my $withinCladeMax = 0.0007126632;
my $withinCladeMax = 0.0002809185;
my $datafile = "dataRopars.LOHnoR.tsv"; # OUTPUT TABLE OF DATA FOR EACH STRAIN FROM HETEROZYGOSITY ANALYSIS
my $hetdata = "backboneRopars.het.tsv";
my $table = "Table2Ropars.tsv";


#my $tiff = "yes";
my $tiff = "no";
my $skip = "yes";
my $colors = "colors4";
my $pew = 2;		# population edge width
my $ew = 2;		# default edge width
my $ec = "#737373"; 	# default edge color

my @backbone = qw(CEC3682 CEC3623 CEC3560 CEC4946 CEC3549 CEC3531 CEC3658 CEC3681 CEC3557 CEC3676 CEC3607 CEC3674 CEC2023 CEC3622 CEC3634 CEC2018 CEC3533 CEC5026 CEC3616 CEC3706 CEC3711 CEC4495 CEC3618 CEC3704 CEC3680 CEC2019 CEC3537 CEC4877 CEC5020 CEC5019 CEC3661 CEC4039 CEC3707 CEC4498 CEC3712 CEC4486 CEC3550 CEC3600 CEC3664 CEC3786 CEC2871 CEC2876);


# includes Clades A, B and 10 diverged strains
my @backbone2 = qw(CEC3682 CEC3623 CEC3560 CEC4946 CEC3549 CEC3531 CEC3637 CEC3681 CEC3557 CEC3676 CEC3607 CEC3674 CEC2023 CEC3622 CEC3634 CEC2018 CEC3533 CEC5026 CEC3616 CEC3706 CEC3711 CEC4495 CEC3618 CEC3704 CEC3680 CEC2019 CEC3537 CEC4877 CEC5020 CEC5019 CEC3661 CEC4039 CEC3707 CEC4498 CEC3712 CEC4486 CEC3550 CEC3600 CEC3664 CEC3786 CEC2871 CEC2876 CEC3548 CEC3715 CEC4038 CEC3708 CEC4254 CEC3555 CEC3613 CEC3619 CEC3663 CEC3665 CEC3671 CEC4261 CEC4481 CEC4485
CEC723);

my %t2oakclades = (
	"NCYC4145" => "unknown",
	"NCYC4144" => "clade 18",
	"NCYC4146" => "clade 4",
	
);

## INFORMATION FOR COLORING PHYLOGENY
my %node = (
	1 => 69,
	2 => 93,
	3 => 114,
	4 => 83,
	8 => 79, 	# CEC2023
	9 => 121,	# CEC2018
	10 => 105,	# CEC3616
	11 => 72,	# CEC4495
	12 => 118,	# CEC2019 
	13 => 99,	# CEC4877
	16 => 101,	# CEC3600
	18 => 109,
	"A" => 112,	# CEC3548
	"B1" => 24,	# CEC3708
	"B2" => 25,	# CEC4254
	"C1" => 21,	# CEC3708
	"C2" => 20,	# CEC4039
	"D1" => 31,	# CEC3707
	"D2" => 32,	# CEC4498
	"E1" => 8,	# CEC3707
	"E2" => 9,	# CEC4486
);


## INFORMATION FOR MAKING TABLE 2

my %mtl = (
	"NCYC4144" => "homa",
	"NCYC4145" => "het",
	"NCYC4146" => "het",
	"NCYC597" => "het",
	"WT_SC5314" => "het",
	"12C" => "homa",
	"L26" => "homa",
	"P37005" => "homa",
	"19F" => "homalpha",
	"P78048" => "homalpha",
	"P37037" => "het",
	"P37039" => "het",
	"P57072" => "homalpha",
	"P78042" => "homalpha",
	"P34048" => "het",
	"P87" => "homa",
	"GC75" => "homalpha",
	"P75016" => "het",
	"P75063" => "het",
	"P94015" => "homa",
	"P60002" => "homa",
	"P75010" => "homalpha",
	"P57055" => "het",
	"P76055" => "het",
	"P76067" => "het",
	"SC5314" => "het",
	"1AA" => "het",
);

my %clades;

open CLADES, "$clades" or die "couldn't open $clades : $!";
while (<CLADES>) {
	if (/^(\S+)\s+(\S+)/m) { $clades{$1} = $2; }
}
close CLADES;

# SET UP SCRIPTS FOR RAXML SNP PHYLOGENY


open SHOUT, ">runSNPphylo" or die "couldn't open runSubPhylo : $!";

printphyloheader("runSNPphylo");

`mkdir raxml$outname`;

print SHOUT"
fa2phylip.pl -i $outname.mfa > $outname.phy
raxmlHPC-PTHREADS-AVX -T 4 -f a -x 12345 -p 12345 -m GTRGAMMA -N $bootstraps -s $outname.phy -n $outname -w $dir/raxml$outname/
";

close SHOUT;

print "\nTo generate $outname.mfa and $foutname from 41467_2018_4787_MOESM6_ESM.txt and chr*.mfa alignements, then run runRopars2al.pl in the same directory as those files\n\n";


print "\nNext step:\nqsub runSNPphylo\n\n";



########################################################
### PAINT snps SHOWING ONLY COLORS with within clade similarity (with -M option: 0.1%)
#######################################################


`mkdir paintWithOakprobeM90`;

open SHOUT, ">runSNPfaChrompaintM90" or die "couldn't open runSNPfaChrompaintM90 : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runSNPfaChrompaintM90
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=64gb

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3

cd  /lustre1/doudab/candidAnalysis/final2/paintWithOakprobeM90
ln -s ../roparsSuppTable3.mfa .

for strain in NCYC4144 NCYC4145 NCYC4146 NCYC597 SC5314
do
  faChrompaint.pl -i roparsSuppTable3.mfa -W 5000 -r \$strain -c ../clades -C ../$colors -e 'NCYC597 1AA SC5314_A' -M $withinCladeMax
done


";

close SHOUT;

print "\nNext step:\nqsub runSNPfaChrompaintM90\n\n";



#######################################################
### Heterozygosity of SNP analysis
#######################################################
my $run = "short";

if ($run eq "long") {


	my $basecomp = `basecomp.pl -i $foutname`;

	print "Ran basecomp.pl .. now reading the output:\n";

	while ($basecomp =~ /(\S+)\s+A : \d+\s+.*?totalACGT : (\d+)\s+.*?totalALL : (\d+)/g) {
		my $strain = $1;
		my $het = ($3-$2)/$3;
		print "$strain\theterozygosity: $het\n";
	}	

}

#######################################################
### SET UP SCRIPTS FOR WHOLE-GENOME MAPPING AND BASECALLING
#######################################################

# 
# Read in description of experiments for each strain name

my (%strains,%readnames);

open IN, "<$ncbireaddesc" or die "couldn't open $ncbireaddesc : $!";

while (<IN>) {
	if (/^(S\S+)\tIllumina.*?strain\s+(\S+),/m) {
		my $expt = $1;
		my $strain = $2;
		$strains{$expt} = $strain;
	#	`mkdir $strain`;
	}
}
close IN;

# Read in description of names of read files for each experiment

open IN, "<$ebireaddesc" or die "couldn't open $ebireaddesc : $!";

while (<IN>) {
	if (/^\S+\s+\S+\s+\S+\s+(\S+).*?ftp.sra.ebi.ac.uk.*\/(\S+)_1.fastq.gz;.*?\S+_2.fastq.gz/m) {
		$readnames{$strains{$1}} = $2;  
		print "strain: $strains{$1}\t$readnames{$strains{$1}}\n";
	#	`ln -s /lustre1/doudab/candidAnalysis/Ropars_etal18/$2*.fastq.gz $strains{$1}/.`;
		
	}
}
close IN;


#open SHOUT, ">runMapCallRopars" or die "couldn't open runMapCallRopars : $!";
open SHOUT, ">runMapCallBatch16" or die "couldn't open runMapCallBatch16 : $!";


print SHOUT "
#PBS -S /bin/bash
#PBS -q batch
#PBS -N mapCallBatch16
#PBS -l nodes=1:ppn=1:dbnode		# reduce nodes because samtools takes ages and only uses 1 processor
#PBS -l walltime=168:00:00
#PBS -l mem=16gb			# job history on these does not suggest very high memory

#PBS -M doudab\@uga.edu 
#PBS -m ae


module load bwa/0.7.10
module load samtools/1.2
module load freebayes/1.0.1
module load R/3.2.3
module load seqtk/1.0-r82-dirty

";


my @batch1 = qw(CEC3555 CEC3556 CEC3557 CEC3558 CEC3559 CEC3560 CEC3561 CEC3579 CEC3596 CEC3597);	# 10 strains
my @batch2 = qw(CEC3600 CEC3601 CEC3602 CEC3603 CEC3605 CEC3607 CEC3609 CEC3610 CEC3611 CEC3612 CEC3613 CEC3614 CEC3615 CEC3616 CEC3617 CEC3618 CEC3619 CEC3621 CEC3622 CEC3623 CEC3626 CEC3627 CEC3634 CEC3636 CEC3637); 	# 25 strains
# CEC3621 CEC3622 CEC3623 CEC3626 CEC3627 CEC3634 CEC3636 CEC3637
my @batch3 = qw(CEC3638 CEC3658 CEC3659 CEC3660 CEC3661 CEC3662 CEC3663 CEC3664 CEC3665 CEC3668 CEC3669 CEC3671 CEC3672 CEC3673 CEC3674 CEC3675 CEC3676 CEC3678 CEC3679 CEC3680);
#  CEC3678 CEC3679 CEC3680
my @batch4 = qw(CEC3681 CEC3682 CEC3683 CEC3685 CEC3686 CEC3704 CEC3706 CEC3707 CEC3708 CEC3711 CEC3712 CEC3715 CEC3716 CEC3786 CEC4032 CEC4035 CEC4038 CEC4104 CCEC4039 CEC4103 EC4106 CEC4107 CEC4108 CEC4254 CEC4256);
# CEC4104 CEC4106 CEC4107 CEC4108 CEC4254 CEC4256
my @batch5 = qw(CEC4259 CEC4261 CEC4479 CEC4480 CEC4481 CEC4482 CEC4484 CEC4485 CEC4486 CEC4487 CEC4488 CEC4489 CEC4492 CEC4493 CEC4494 CEC4495 CEC4496 CEC4497 CEC4498 CEC4499 CEC4500 CEC4525 CEC4854 CEC4855 CEC4856); 
my @batch6 =qw(CEC4857 CEC4858 CEC4859 CEC4860 CEC4861 CEC4864 CEC4865 CEC4877 CEC4878 CEC4879 CEC4881 CEC4882 CEC4883 CEC4884 CEC4886 CEC4887 CEC4888 CEC4889 CEC4943 CEC4944 CEC4945 CEC4946);
#CEC4944 CEC4945 CEC4946)
my @batch7 =qw(CEC5019 CEC5020 CEC5021 CEC5022 CEC5023 CEC5024 CEC5025 CEC5026 CEC5027 CEC5028 CEC5029 CEC5030 CEC704 CEC708 CEC709 CEC711 CEC712 CEC716 CEC718 CEC723);
# CEC718 CEC723

# second round:
my @batch8 = qw(CEC3553 CEC3658 CEC3621 CEC3622 CEC3623 CEC3626 CEC3627 CEC3634 CEC3636 CEC3637 CEC3678 CEC3679 CEC3680 CEC4104 CCEC4039 CEC4103 CEC4106);
my @batch9 = qw(CEC4855 CEC4856 CEC4107 CEC4108 CEC4254 CEC4256 CEC4497 CEC4498 CEC4499 CEC4500 CEC4525 CEC4854 CEC718 CEC723);

# third round
my @batch10 = qw(CEC4884 CEC4886 CEC4887 CEC4888);
my @batch11 = qw(CEC4889 CEC4943 CEC4944 CEC4945);
my @batch12 = qw(CEC4946 CEC718 CEC723 CEC711); 
my @batch13 = qw(CEC712 CEC716);

my @batch14 = qw(CEC4495);
my @batch15 = qw(CEC4498);
my @batch16 = qw(CEC4946);


#foreach my $strain (sort values %strains) {
foreach my $strain (@batch16) {

	print SHOUT "

cd /lustre1/doudab/candidAnalysis/final2/$strain

bwa mem -t 1 /lustre1/doudab/candidAnalysis/final/SC5314_A22edited.mfa $readnames{$strain}_1.fastq.gz $readnames{$strain}_2.fastq.gz | samtools view -b - | samtools sort - $strain

samtools mpileup -I -d 10000 -uf /lustre1/doudab/candidAnalysis/final/SC5314_A22edited.mfa $strain.bam | bcftools call -c > $strain.I.vcf

# create a link for subsequent vcf analysis in the main directory
ln -s /lustre1/doudab/candidAnalysis/final2/$strain/$strain.I.vcf /lustre1/doudab/candidAnalysis/final2/$strain.I.vcf

# VISUALIZE THE VCF FILE						# This must not have worked (no gff file) in final2
ln -s /lustre1/doudab/candidAnalysis/final/SC5314annotations.gff .
vcf2allelePlot.pl -i $strain.I.vcf -g SC5314annotations.gff -h -m

# CONVERT VCF TO FASTA
vcfutils.pl vcf2fq $strain.I.vcf > $strain.I.fq
 seqtk seq -q 40 -A $strain.I.fq >$strain.I.fa

# CREATE ALIGNMENTS FOR EACH CHROMOSOME
for chr in chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chrR
do
  faChoose.pl -i $strain.I.fa -n \$chr -o $strain\$chr.I.fa 
  cat $strain\$chr.I.fa | sed s/\$chr/$strain/ > temp
  mv temp $strain\$chr.I.fa
done
echo \"created 8 $strain chr*.I.fa files\"
";


}

#print "\nNext step:\nqsub runMapCallRopars\n";
print "\nNext step:\nqsub runMapCallBatch16\n";


close SHOUT;	

#######################################################
# CREATE CHROMOSOME ALIGNMENTS
#######################################################


open SHOUT, ">runChrAlign" or die "couldn't open runChrAlign : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runChrAlign
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=168:00:00
#PBS -l mem=64gb

#PBS -M doudab\@uga.edu 
#PBS -m ae


cd  /lustre1/doudab/candidAnalysis/final2/

";


# Prep whole chromosome alignment files 

my ($allchrfiles,$allname);

foreach my $chr (@chr) {



	$allchrfiles .= "$chr.mfa "; $allname .= $chr;

	print SHOUT "
cat ../final/NCYC*/*NCYC*$chr.I.fa > $chr.mfa
cat */*$chr.I.fa >> $chr.mfa
fastaLC2n.pl -i $chr.mfa -o temp
alcat.pl -i temp -o $chr.mfa
fa2phylip.pl -i $chr.mfa > $chr.phy

";
	
}

close SHOUT;

print "\nNext step:\nqsub runChrAlign\n\n";


#######################################################
### SET UP SCRIPTS FOR chromosome painting ALL STRAINS (no backbone) without divergence filter
#######################################################



# Painting


`mkdir paintNoNCYCprobe`;
`ln -s /lustre1/doudab/candidAnalysis/final2/chr*.mfa /lustre1/doudab/candidAnalysis/final2/paintNoNCYCprobe/.`;

my @pbatch1 = qw(CEC3555 CEC3556 CEC3557 CEC3558 CEC3559 CEC3560 CEC3561 CEC3579 CEC3596 CEC3597);	# 10 strains
my @pbatch2 = qw(CEC1289 CEC1424 CEC1426 CEC1427 CEC1492 CEC2018 CEC2019 CEC2020 CEC2021 CEC2022);
my @pbatch3 = qw(CEC2023 CEC2871 CEC2872 CEC2876 CEC3494 CEC3529 CEC3530 CEC3531 CEC3532 CEC3533);
my @pbatch4 = qw(CEC3534 CEC3535 CEC3536 CEC3537 CEC3540 CEC3541 CEC3544 CEC3547 CEC3548 CEC3549);
my @pbatch5 = qw(CEC3550 CEC3551 CEC3553 CEC3658 CEC3600 CEC3601 CEC3602 CEC3603 CEC3605 CEC3607);
my @pbatch6 = qw(CEC3609 CEC3610 CEC3611 CEC3612 CEC3613 CEC3614 CEC3615 CEC3616 CEC3617 CEC3618);
my @pbatch7 = qw(CEC3619 CEC3621 CEC3622 CEC3623 CEC3626 CEC3627 CEC3634 CEC3636 CEC3637 CEC3638);
my @pbatch8 = qw(CEC3658 CEC3659 CEC3660 CEC3661 CEC3662 CEC3663 CEC3664 CEC3665 CEC3668 CEC3669);

open SHOUT, ">runpaintNCYC" or die "couldn't open runpaintNCYC : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runpaintNCYC
#PBS -l nodes=1:ppn=1:dbnode 
#PBS -l walltime=168:00:00
#PBS -l mem=16gb

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3

cd  /lustre1/doudab/candidAnalysis/final2/paintNoNCYCprobe

for strain in NCYC4144 NCYC4145 NCYC4146 NCYC597 

#for strain in @pbatch2

do
  faChrompaint.pl -I chr -r \$strain -c ../clades -C ../$colors -e NCYC
done


";

close SHOUT;

print "\nNext step:\nqsub runpaintNCYC\n\n";



#######################################################
### SET UP SCRIPTS FOR chromosome painting ONLY USING BACKBONE without divergence filter
#######################################################

open SHOUT, ">runBChrAlign" or die "couldn't open runBChrAlign : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runBChrAlign
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=168:00:00
#PBS -l mem=64gb

#PBS -M doudab\@uga.edu 
#PBS -m ae


cd  /lustre1/doudab/candidAnalysis/final2/

";



# Prep whole chromosome alignment files 

$allchrfiles = "";
$allname = "";

print "Preparing script for making alignments for NCYC strains, SC5314 and ".@backbone." backbone strains\n"; 

foreach my $chr (@chr) {

	$allchrfiles .= "b$chr.mfa "; $allname .= $chr;

	print SHOUT "
cat ../final/NCYC*/*NCYC*$chr.I.fa > b$chr.mfa
cat ../final/SC5314/SC5314$chr.I.fa >> b$chr.mfa
";

	foreach my $bstrain (@backbone2) { print SHOUT "cat $bstrain/$bstrain$chr.I.fa >> b$chr.mfa\n"; }

	print SHOUT "
fastaLC2n.pl -i b$chr.mfa -o temp
alcat.pl -i temp -o b$chr.mfa
fa2phylip.pl -i $chr.mfa > b$chr.phy
";
	
}

print SHOUT "
alcat.pl -i '$allchrfiles' -o $allname.mfa
fa2phylip.pl -i $allname.mfa > $allname.phy
";

close SHOUT;

print "\nNext step:\nqsub runBChrAlign\n\n";


# Painting
# change to @backbone2

#my @backbone2 = qw(CEC3682 CEC3623 CEC3560 CEC4946 CEC3549 CEC3531 CEC3637 CEC3681 CEC3557 CEC3676 CEC3607 CEC3674 CEC2023 CEC3622 CEC3634 CEC2018 CEC3533 CEC5026 CEC3616 CEC3706 CEC3711 CEC4495 CEC3618 CEC3704 CEC3680 CEC2019 CEC3537 CEC4877 CEC5020 CEC5019 CEC3661 CEC4039 CEC3707 CEC4498 CEC3712 CEC4486 CEC3550 CEC3600 CEC3664 CEC3786 CEC2871 CEC2876 CEC3548 CEC3715 CEC4038 CEC3708 CEC4254 CEC3555 CEC3613 CEC3619 CEC3663 CEC3665 CEC3671 CEC4261 CEC4481 CEC4485
#CEC723);


#`mkdir paintBNoNCYCprobe`;
`mkdir paintB2NoNCYCprobe`;

`ln -s /lustre1/doudab/candidAnalysis/final2/bchr*.mfa /lustre1/doudab/candidAnalysis/final2/paintB2NoNCYCprobe/.`;


my @bbatch1 = qw(CEC2019 CEC3537 CEC4877 CEC5020 CEC5019 CEC3661 CEC4039 CEC3707);
my @bbatch2 = qw(CEC4498 CEC3712 CEC4486 CEC3550 CEC3600 CEC3664 CEC3786 CEC2871);
my @bbatch3 = qw(CEC2876 CEC3548 CEC3715 CEC4038 CEC3708 CEC4254 CEC3555 CEC3613);
my @bbatch4 = qw(CEC3619 CEC3663 CEC3665 CEC3671 CEC4261 CEC4481 CEC4485);
my @bbatch5 = qw(CEC3682 CEC3623 CEC3560 CEC4946 CEC3549 CEC3531 CEC3637 CEC3704 CEC3680);
my @bbatch6 = qw(CEC3681 CEC3557 CEC3676 CEC3607 CEC3674 CEC2023 CEC3622 CEC3634);
my @bbatch7 = qw(CEC2018 CEC3533 CEC5026 CEC3616 CEC3706 CEC3711 CEC4495 CEC3618);



open SHOUT, ">runpaintB2Backbone7" or die "couldn't open runpaintB2Backbone7 : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runpaintB2Batch7
#PBS -l nodes=1:ppn=1:dbnode 
#PBS -l walltime=168:00:00
#PBS -l mem=16gb

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3

cd  /lustre1/doudab/candidAnalysis/final2/paintB2NoNCYCprobe

#for strain in @backbone2 

for strain in @bbatch7

do
  faChrompaint.pl -I bchr -r \$strain -c ../clades -C ../$colors -e 'NCYC SC5314'
done


";

close SHOUT;

print "\nNext step:\nqsub runpaintB2Backbone7\n\n";


## FIGURE OUT 90% threshold for within clade divergence and visualise together with between clade divergence

open RCMD, ">paintB2NoNCYCprobe/divHist.R" or die "couldn't open paintB2NoNCYCprobe/divHist.R : $!";

print RCMD "

pdf('divHist.pdf')
options(scipen=999)

par(mfrow=c(2,1))
sameclade<-read.table(\"sameclade.pdiffs\")
head(sameclade)

quantile(sameclade\$V9,c(0.9,0.95,0.99),na.rm=T)
summary(sameclade\$V9)

hist(sameclade\$V9,0:170/10000,col=\"blue\",xlab=\"Proportion of bases differing between pairs of sequences\",main=\"Within-clade comparisons\")
#abline(v=quantile(sameclade\$V9,0.95,na.rm=T),col=\"red\")
abline(v=quantile(sameclade\$V9,0.9,na.rm=T),col=\"green\")

diffclade<-read.table(\"diffclade.pdiffs\")
summary(diffclade\$V9)
head(diffclade)
quantile(diffclade\$V9,c(0.05,0.1,0.15,0.2),na.rm=T)

hist(diffclade\$V9,0:170/10000,col=\"purple\",xlab=\"Proportion of bases differing between pairs of sequences\",main=\"Between-clade comparisons\")
abline(v=quantile(sameclade\$V9,0.9,na.rm=T),col=\"green\")
#abline(v=quantile(sameclade\$V9,0.95,na.rm=T),col=\"red\")
";
close RCMD;

open SHOUT, ">runDivhist" or die "couldn't open runDivhist : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runDivhist
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=64gb

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3

cd  /lustre1/doudab/candidAnalysis/final2/paintB2NoNCYCprobe

cat *.pDiffs.tsv | grep -v NC | grep chr | awk '\$3==\$5{ print \$0; }' > sameclade.pdiffs  # only within clade diffs because \$3==\$5: oak strains and singletons excluded
cat *.pDiffs.tsv | grep -v NC | grep -v pDiff | grep chr | awk '\$5!=\$3 { print \$0; }' > diffclade.pdiffs

`R CMD BATCH --no-save --no-restore divHist.R`
";


close SHOUT;

print "\nNext step:\nqsub runDivhist\n\n";


########################################################
### PAINT CHROMOSOMES SHOWING ONLY COLORS with within clade similarity (<90% with -M option)
#######################################################


`mkdir paintWithOakprobeM00066B`;
`ln -s /lustre1/doudab/candidAnalysis/final2/bchr*.mfa /lustre1/doudab/candidAnalysis/final2/paintWithOakprobeM00066B/.`;


open SHOUT, ">runfaChrompaintM00066B" or die "couldn't open runfaChrompaintM00066B : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runfaChrompaintM90
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=168:00:00
#PBS -l mem=16gb

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3

cd  /lustre1/doudab/candidAnalysis/final2/paintWithOakprobeM00066B

#for strain in NCYC4144 NCYC4145 NCYC4146 NCYC597 SC5314 @backbone
for strain in NCYC4144 NCYC4145 NCYC4146 NCYC597 SC5314 @backbone2
do
  faChrompaint.pl -I bchr -r \$strain -c ../clades -C ../$colors -e 'NCYC597 SC5314' -M $withinCladeMax
done


";

close SHOUT;

print "\nNext step:\nqsub runfaChrompaintM00066B\n\n";


########################################################
### PAINT Genome SHOWING ONLY COLORS with within clade similarity (with -M option: 0.1%)
#######################################################


`mkdir paintWithOakGenomeM90`;
`ln -s /lustre1/doudab/candidAnalysis/final2/bchr*.mfa /lustre1/doudab/candidAnalysis/final2/paintWithOakGenomeM90/.`;


open SHOUT, ">runGenomefaChrompaintM90" or die "couldn't open runGenomefaChrompaintM90 : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runGenomefaChrompaintM90
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=168:00:00
#PBS -l mem=16gb

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3

cd  /lustre1/doudab/candidAnalysis/final2/paintWithOakGenomeM90

for strain in NCYC4144 NCYC4145 NCYC4146 NCYC597 SC5314 @backbone2
do
  faChrompaint.pl -I bchr -r \$strain -c ../clades -C ../$colors -e 'NCYC597 SC5314' -M $withinCladeMax
done


";

close SHOUT;

print "\nNext step:\nqsub runGenomefaChrompaintM90\n\n";



# SET UP SCRIPTS FOR RAXML GENOME PHYLOGENY


open SHOUT, ">runGenomePhylo2" or die "couldn't open runGenomePhylo2 : $!";

printphyloheader("runGenomePhylo2");

$outname = "genome2";

`mkdir raxml$outname`;

print SHOUT"
raxmlHPC-PTHREADS-AVX -T 16 -f a -x 12345 -p 12345 -m GTRGAMMA -N $bootstraps -s chr1chr2chr3chr4chr5chr6chr7chrR.phy -n $outname -w $dir/raxml$outname/
";

close SHOUT;


print "# SET UP SCRIPTS FOR RAXML GENOME PHYLOGENY\nNext step:\nqsub runGenomePhylo2\n\n";



## DRAW 1-PART PHYLOGENETIC FIGURE : genome-wide

$outname = "genome2";

my %colors;

open COLORS, "$colors" or die "couldn't open $colors : $!";
while (<COLORS>) {
	if (/^(\S+)\s+(\S+)/m) { 
		$colors{$1} = $2; 
	}
}

open RCMD, ">phyloRopars.R" or die "couldn't open phyloRopars.R : $!";

print RCMD "library(ape)\n";

if ($tiff eq "no") { print RCMD "pdf('phyloRopars.pdf',7,9)\n"; }
else { print RCMD "bitmap(\"phyloRopars.tiff\", type = \"tiff24nc\", res = 300)\n"; }

print RCMD "

par(mfrow=c(1,1),mar=c(3,3,3,3),cex=1,xpd=T)
";

printRcmds4fig($outname);

print RCMD "
par(mar=c(0,1,0,1),cex=1)
#legend(\"bottomleft\", \"a.\", bty=\"n\",cex=2)
#legend(\"bottomleft\",c(\"oak\",\"clade 1\",\"clade 2\",\"clade 3\",\"clade 4\",\"clade 8\",\"clade 9\",\"clade 10\",\"clade 11\",\"clade 12\",\"clade 13\",\"clade 16\",\"clade 18\",\"clade A\",\"clade B\",\"clade C\",\"clade D\",\"clade E\",\"unknown\"),text.col=c(wildcol,clade1col,clade2col,clade3col,clade4col,clade8col,clade9col,clade10col,clade11col,clade12col,clade13col,clade16col,clade18col,cladeAcol,cladeBcol,cladeCcol,cladeDcol,cladeEcol,unknowncol),bty=\"n\",cex=1.1)

legend(\"bottomleft\",c(\"oak\",\"unknown\"),text.col=c(wildcol,unknowncol),bty=\"n\",cex=1.1)

#,inset = -0.05

text(0.008,9,\"clade 1\",col=clade1col)
text(0.00925,33,\"clade 2\",col=clade2col)
text(0.00805,47,\"clade 3\",col=clade3col)
text(0.0089,21.6,\"clade 4\",col=clade4col)
text(0.0073,18,\"clade 8\",col=clade8col)
text(0.0078,3,\"clade 9\",col=clade9col)
text(0.0072,37,\"clade 10\",col=clade10col)
text(0.0087,12,\"clade 11\",col=clade11col)
text(0.0073,51,\"clade 12\",col=clade12col)
text(0.0068,58,\"clade 13\",col=clade13col)
text(0.0041,61,\"clade 16\",col=clade16col)
text(0.00805,44,\"clade A\",col=cladeAcol)
text(0.0088,30.5,\"clade B\",col=cladeBcol)
text(0.0087,27,\"clade C\",col=cladeCcol)
text(0.0047,55.5,\"clade D\",col=cladeDcol)
text(0.0072,14.5,\"clade E\",col=cladeEcol)

text(0.00805,41,\"clade 18\",col=clade18col)

";

close RCMD;

print "\nPrinted R commands for drawing 1-part Figure that shows genome-wide phylogeny only\n";

open SHOUT, ">runPhyloFigRopars" or die "couldn't open runPhyloFigRopars : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runPhyloFigRopars
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=16gb

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3

cd /lustre1/doudab/candidAnalysis/final2/

`R CMD BATCH --no-save --no-restore phyloRopars.R`

";


if ($tiff eq "yes") { print SHOUT "`tiff2pdf phyloRopars.tiff -o phyloRopars.pdf`\n"; }

close SHOUT;

print "\nNext step to run the R commands:\nqsub runPhyloFigRopars\n\n";




## CREATE LOH FILTER AND MORE VCF2ALLELES PLOTS

unless ($skip eq "yes") {
	`ln -s /lustre1/doudab/candidAnalysis/final/SC5314annotations.gff /lustre1/doudab/candidAnalysis/final2/.`;
	`ln -s /lustre1/doudab/candidAnalysis/final/NCYC*.vcf /lustre1/doudab/candidAnalysis/final2/.`;
	`ln -s */*.I.vcf .`; 
}

opendir (DIR, ".") or die "couldn't open . : $!";

my (@vcffiles,%simple2dhetfiles);						# collect the names of data directories, strains and their study accessions
while (defined (my $name = readdir(DIR))) {

	if ($name =~ /^(\S+)\.I\.vcf$/i) { 
#		print "$name\tstrain:$1\n"; 	
		$strains{$name} = $1;
		push (@vcffiles, $name);
	}

	if ($name =~ /^(\S+)\.simple2d\.het\.txt$/i) { 		# already analysed?
		$simple2dhetfiles{"$1.I.vcf"}++;
	}
}		



# estimate heterozygosity for every strain after filtering out LOH regions
print "# The following scripts will visualize the vcf plots and\n";
print "# estimate heterozygosity for every strain after filtering out LOH regions\n";


print "\n# Next step:\n\n";


foreach my $file (sort @vcffiles) { 

	if (exists $simple2dhetfiles{$file}) { next; }

	open SHOUT, ">runVcfPlots$strains{$file}" or die "couldn't open runVcfPlots$strains{$file} : $!";
	print SHOUT "
#PBS -S /bin/bash
#PBS -q batch
#PBS -N vcf2Plots$strains{$file}
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=8gb		# these jobs usually only require 5GB of memory, and increasing memory increases wait in the queue

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3		# vcf2allelePlot.pl uses R



cd /lustre1/doudab/candidAnalysis/final2

# VISUALIZE THE VCF FILES
";

	print SHOUT "vcf2allelePlot.pl -i $strains{$file}.I.vcf -g SC5314annotations.gff -D\n"; 



	print SHOUT "cp $strains{$file}.I.gff $strains{$file}.LOH.I.gff\n"; 
	
	# filter chrR from heterozygosity estimates
	print SHOUT "echo \"chrR\tvcf2alleleplot.pl\tchrR\t0\t2300000\t.\t+\t.\tfilter chrR which differs in ploidy\" >> $strains{$file}.LOH.I.gff\n"; 

	print SHOUT "\n# FILTER LOH REGIONS\n\n"; 			# USING LOH ANNOTATIONS CREATED ABOVE (REMOVE COMMENTS above to regenerate LOH annotations)


	print SHOUT "cp $strains{$file}.I.pdf $strains{$file}.simplenoLOH.pdf\n"; 
	print SHOUT "vcf2allelePlot.pl -i $strains{$file}.I.vcf -g $strains{$file}.LOH.I.gff -o $strains{$file}.simple2d -D\n"; 

	close SHOUT;	

	print "qsub runVcfPlots$strains{$file}\n";

}

###########################################################
## CODE ADAPTED FROM runHetanalysis.pl for TABLE 2 AND HETEROZYGOSITY ANALYSES
###########################################################

open TABLE, ">$table" or die "couldn't open $table: $!";

opendir (DIR, ".") or die "couldn't open . : $!";

#my (%strains,@hetfiles);					# collect the names of het summary files and strains 
my (@hetfiles);					# collect the names of het summary files and strains 

`ln -s ../final/SC5314.simple2d.het.txt .`;	# include SC5314 in the heterozygosity analysis
`ln -s ../final/SC5314.LOH.I.gff .`;		# include SC5314 in the heterozygosity analysis

while (defined (my $name = readdir(DIR))) {

	if ($name =~ /^(\S+)\.simple2d\.het\.txt$/i) { 
	#	print "$name\tstrain:$1\n"; 	
		$strains{$name} = $1;
		push (@hetfiles, $name);
	}

}		

# EXTRACT GENOMEWIDE HETEROZYGOSITY DATA FROM *\.simple2d\.het\.txt files


open DATA, ">$datafile" or die "couldn't open $datafile : $!";
print DATA "strain\tfile\tclade\tq40ps\tq40l\tq40het\tfps\tfl\tfhet\tfdps\tfdl\tfdhet\tlohl\n";
print TABLE "Strain or Clade & MTL & Mean heterozygosity \\newline (max. - min.)\\footnote\\ & Filtered length (Mbp)\\footnote\\ & Mean filtered heterozygosity\\footnote\\ \\newline (max. - min.) \\\\ \\hline\n";


my (%hqh2strains,%hqh,%lohl,%dfl,%dfh,%dfps,%hqps,%hql,%fl,%fps,%fh);

my $nncyc = 0;
my $sumhqhncyc = 0;
my $sumhqhclinic = 0;
my $nclinic = 0;
my $sumlohncyc = 0;
my $sumlohclinic = 0;
my $sumdfhncyc = 0;
my $sumdfhclinic = 0;
foreach my $file (sort @hetfiles) { 

	unless (defined $mtl{$strains{$file}}) { $mtl{$strains{$file}} = "?"; }
	unless (exists $clades{$strains{$file}}) { $clades{$strains{$file}} = "?"; }

	print DATA "$strains{$file}\t$file\t$clades{$strains{$file}}";

	open IN, "$file" or die "couldn't open $file : $!";
	my ($fl,$fh);
	while (<IN>) {


		if (/(\S+)\s+\(\s+(\d+)\s+\/\s+(\d+)\s+\)\s+# Genomewide heterozygosity/) {
		#	print "genomewide $strains{$file}:\t$1\t$2\t$3\n";
			print DATA "\t$2\t$3\t$1";
			$hqh{$strains{$file}} = $1;
			$hqps{$strains{$file}} = $2;
			$hql{$strains{$file}} = $3;

			$hqh2strains{$1} = $strains{$file}; 	# have a hash where the keys can be sorted by q40het
								# this will only work because there are no ties
			if ($strains{$file} =~ /NCYC4/) { $sumhqhncyc += $hqh{$strains{$file}}; $nncyc++; }
			elsif ($strains{$file} eq "1AA") { print "skipping 1AA from estimate of genome-wide means\n"; }
			else { $sumhqhclinic += $hqh{$strains{$file}}; $nclinic++; }
		}
		if (/(\S+)\s+\(\s+(\d+)\s+\/\s+(\d+)\s+\)\s+# Unannotated heterozygosity/) {
		#	print "unannotated $strains{$file}:\t$1\t$2\t$3\n";
			print DATA "\t$2\t$3\t$1";
			$fh{$strains{$file}} = $1;
			$fps{$strains{$file}} = $2;
			$fl{$strains{$file}} = $3;

		}
		if (/(\S+)\s+\(\s+(\d+)\s+\/\s+(\d+)\s+\)\s+# Depth <=/) {
		#	print "$strains{$file}:\t$1\t$2\t$3\n";
			print DATA "\t$2\t$3\t$1";
			$dfh{$strains{$file}} = $1;
			$dfps{$strains{$file}} = $2;
			$dfl{$strains{$file}} = $3;

			if ($strains{$file} =~ /NCYC4/) { $sumdfhncyc += $dfh{$strains{$file}}; }
			elsif ($strains{$file} eq "1AA") { print "skipping 1AA from estimate of genome-wide means\n"; }
			else { $sumdfhclinic += $dfh{$strains{$file}}; }
		}
	}	
	close IN;

	my $gffile = "$strains{$file}.LOH.I.gff";
	open LOH, "<$gffile" or die "couldn't open $gffile : $!";
	my $loh;
	while (<LOH>) {
		if (/chr\S+\s+vcf2allelePlot\.pl\s+LOH\s+(\d+)\s+(\d+)/) {
			$loh += $2-$1;
		#	print "$strains{$file}\t$1\t$2\t$loh\n";
		}		
	}
#	print "extracted LOH data for $strains{$file} from $gffile.\n";
	if (defined $loh) { 
		print DATA "\t$loh\n";
		$lohl{$strains{$file}} = $loh;
		if ($strains{$file} =~ /NCYC4/) { $sumlohncyc += $loh; }
		elsif ($strains{$file} eq "1AA") { print "skipping 1AA from estimate of genome-wide means\n"; }
		else { $sumlohclinic += $loh; } 
	}
	else { 
		print "No LOH for $strains{$file} in $gffile.\n"; 
		$lohl{$strains{$file}} = 0; 
		print DATA "\t0\n"; 
	}
	close LOH;
}

close DATA;

print "\nPrinted a summary of heterozygosity data for all ".@hetfiles." *.simple2d.het.txt files to $datafile.\n";




#####################################################################
# EXTRACTING THE DATA FOR TABLE 2 FROM Hetanalysis.Rcmd.Rout 
#####################################################################

open HETANALYSIS, "<Hetanalysis.Rcmd.Rout" or die "couldn't open Hetanalysis.Rcmd.Rout : $!";

my $hetanalysis;
while (<HETANALYSIS>) { $hetanalysis .= $_; }
close HETANALYSIS;

#print $hetanalysis;

my ($meanoakfdl,$meanclinicfdl,$mincfdhet,$maxcfdhet,$minchet,$maxchet,$minc4het,$maxc4het,$meanc4het,$minc18het,$maxc18het,$meanc18het,$mincNChet,$maxcNChet,$meancNChet,$ncc4,$ncc18,$nccu,$minc4fdl,$maxc4fdl,$meanc4fdl,$minc18fdl,$maxc18fdl,$meanc18fdl,$mincNCfdl,$maxcNCfdl,$meancNCfdl,$minc4fdhet,$maxc4fdhet,$meanc4fdhet,$minc18fdhet,$maxc18fdhet,$meanc18fdhet,$mincNCfdhet,$maxcNCfdhet,$meancNCfdhet);

if ($hetanalysis =~ /summary\(oakfdl.*?\d+\s+\d+\s+\d+\s+(\d+)/s) { $meanoakfdl = $1/1000000; }
if ($hetanalysis =~ /summary\(clinicfdl.*?\d+\s+\d+\s+\d+\s+(\d+)/s) { $meanclinicfdl = $1/1000000; }
if ($hetanalysis =~ /summary\(clinicfdhet.*?\s+([0-9\.]+)\s+[0-9\.]+\s+[0-9\.]+\s+[0-9\.]+\s+[0-9\.]+\s+([0-9\.]+)/s) { $mincfdhet = $1; $maxcfdhet = $2; }
if ($hetanalysis =~ /summary\(clinicq40het.*?\s+([0-9\.]+)\s+[0-9\.]+\s+[0-9\.]+\s+[0-9\.]+\s+[0-9\.]+\s+([0-9\.]+)/s) { $minchet = $1; $maxchet = $2; }
if ($hetanalysis =~ /summary\(q40het\[clade==4&strain!="NCYC4146"\].*?\s+([0-9\.]+)\s+[0-9\.]+\s+[0-9\.]+\s+([0-9\.]+)\s+[0-9\.]+\s+([0-9\.]+)/s) { $minc4het = $1; $meanc4het = $2; $maxc4het = $3; }
if ($hetanalysis =~ /summary\(q40het\[clade==18&strain!="NCYC4144"\].*?\s+([0-9\.]+)\s+[0-9\.]+\s+[0-9\.]+\s+([0-9\.]+)\s+[0-9\.]+\s+([0-9\.]+)/s) { $minc18het = $1; $meanc18het = $2; $maxc18het = $3; }
if ($hetanalysis =~ /summary\(q40het\[clade=="NC"&strain!="NCYC4145"\].*?\s+([0-9\.]+)\s+[0-9\.]+\s+[0-9\.]+\s+([0-9\.]+)\s+[0-9\.]+\s+([0-9\.]+)/s) { $mincNChet = $1; $meancNChet = $2; $maxcNChet = $3; }
if ($hetanalysis =~ /length\(q40het\[clade==4&strain!="NCYC4146"\].*?\[1\]\s+(\d+)/s) { $ncc4 = $1; }
if ($hetanalysis =~ /length\(q40het\[clade==18&strain!="NCYC4144"\].*?\[1\]\s+(\d+)/s) { $ncc18 = $1; }
if ($hetanalysis =~ /length\(q40het\[clade==\"NC\"&strain!="NCYC4145"\].*?\[1\]\s+(\d+)/s) { $nccu = $1; }

													# extracted min and max fdl even though I don't use them
if ($hetanalysis =~ /summary\(fdl\[clade==4&strain!="NCYC4146"\].*?\s+([0-9\.]+)\s+[0-9\.]+\s+[0-9\.]+\s+([0-9\.]+)\s+[0-9\.]+\s+([0-9\.]+)/s) { $minc4fdl = $1; $meanc4fdl = $2/1000000; $maxc4fdl = $3; }
if ($hetanalysis =~ /summary\(fdl\[clade==18&strain!="NCYC4144"\].*?\s+([0-9\.]+)\s+[0-9\.]+\s+[0-9\.]+\s+([0-9\.]+)\s+[0-9\.]+\s+([0-9\.]+)/s) { $minc18fdl = $1; $meanc18fdl = $2/1000000; $maxc18fdl = $3; }
if ($hetanalysis =~ /summary\(fdl\[clade=="NC"&strain!="NCYC4145"\].*?\s+([0-9\.]+)\s+[0-9\.]+\s+[0-9\.]+\s+([0-9\.]+)\s+[0-9\.]+\s+([0-9\.]+)/s) { $mincNCfdl = $1; $meancNCfdl = $2/1000000; $maxcNCfdl = $3; }

if ($hetanalysis =~ /summary\(fdhet\[clade==4&strain!="NCYC4146"\].*?\s+([0-9\.]+)\s+[0-9\.]+\s+[0-9\.]+\s+([0-9\.]+)\s+[0-9\.]+\s+([0-9\.]+)/s) { $minc4fdhet = $1; $meanc4fdhet = $2; $maxc4fdhet = $3; }
if ($hetanalysis =~ /summary\(fdhet\[clade==18&strain!="NCYC4144"\].*?\s+([0-9\.]+)\s+[0-9\.]+\s+[0-9\.]+\s+([0-9\.]+)\s+[0-9\.]+\s+([0-9\.]+)/s) { $minc18fdhet = $1; $meanc18fdhet = $2; $maxc18fdhet = $3; }
if ($hetanalysis =~ /summary\(fdhet\[clade=="NC"&strain!="NCYC4145"\].*?\s+([0-9\.]+)\s+[0-9\.]+\s+[0-9\.]+\s+([0-9\.]+)\s+[0-9\.]+\s+([0-9\.]+)/s) { $mincNCfdhet = $1; $meancNCfdhet = $2; $maxcNCfdhet = $3; }



#print "hello $meanclinicfdl\n"; 

my $nstrains = scalar keys %hqh2strains;


foreach my $strain qw(NCYC4146 NCYC4144 NCYC4145) {

	my $mtl;
	if ($mtl{$strain} eq "het") { $mtl = "\$a/\\alpha\$"; }
	elsif ($mtl{$strain} eq "homa") { $mtl = "\$a/a\$"; }
	elsif ($mtl{$strain} eq "homalpha") { $mtl = "\$\\alpha/\\alpha\$"; }
	else { $mtl = $mtl{$strain}; }

	print TABLE "$strain, $t2oakclades{$strain} & $mtl & ".sprintf("%.4f",$hqh{$strain})." & ".sprintf("%.1f",$dfl{$strain}/1000000)." & ".sprintf("%.4f",$dfh{$strain})." \\\\ \n";

	if ($strain eq "NCYC4145") { 		# is this the last entry?

		print TABLE "\n\\textbf{3 oak strains} &  & \\textbf{".sprintf("%.4f",$sumhqhncyc/$nncyc)."} & ".sprintf("%.1f",$meanoakfdl)." & \\textbf{".sprintf("%.4f",$sumdfhncyc/$nncyc)."} \\\\ \\\\ \n";
	
		print TABLE "
\\textbf{".($nstrains-6)." clinical strains} from 17 clades  & 176 $a/\alpha$, \newline 3 $\alpha/\alpha$, \newline 1 $a/a$ & \\textbf{".sprintf("%.4f",$sumhqhclinic/$nclinic)."}\\footnote\\ \\newline (".sprintf("%.4f",$minchet)." - ".sprintf("%.4f",$maxchet).") & ".sprintf("%.1f",$meanclinicfdl)." & \\textbf{".sprintf("%.4f",$sumdfhclinic/$nclinic)."}\\footnote\\ \\newline (".sprintf("%.4f",$mincfdhet)." - ".sprintf("%.4f",$maxcfdhet).") \\\\  \\\\ 

$ncc4 clade 4 strains & & ".sprintf("%.4f",$meanc4het)." \\newline (".sprintf("%.4f",$minc4het)." - ".sprintf("%.4f",$maxc4het).") & ".sprintf("%.1f",$meanc4fdl)." & ".sprintf("%.4f",$meanc4fdhet)." \\newline (".sprintf("%.4f",$minc4fdhet)." - ".sprintf("%.4f",$maxc4fdhet).") \\\\ 

$ncc18 clade 18 strains & & ".sprintf("%.4f",$meanc18het)." \\newline (".sprintf("%.4f",$minc18het)." - ".sprintf("%.4f",$maxc18het).") & ".sprintf("%.1f",$meanc18fdl)." & ".sprintf("%.4f",$meanc18fdhet)." \\newline (".sprintf("%.4f",$minc18fdhet)." - ".sprintf("%.4f",$maxc18fdhet).") \\\\ 

$nccu unknown strains & & ".sprintf("%.4f",$meancNChet)." \\newline (".sprintf("%.4f",$mincNChet)." - ".sprintf("%.4f",$maxcNChet).") & ".sprintf("%.1f",$meancNCfdl)." & ".sprintf("%.4f",$meancNCfdhet)." \\newline (".sprintf("%.4f",$mincNCfdhet)." - ".sprintf("%.4f",$maxcNCfdhet).") \\\\ \\hline
";

	}


}
close TABLE;

print "\nPrinted data for all $nstrains strains from ".@hetfiles." *.simple2d.het.txt files to $table.\nUsed information from Hetanalysis.Rcmd.Rout to make $table. Make sure Hetanalysis.Rcmd.Rout is up-to-date!\n\n";


my %backbone2;
foreach my $strain (@backbone2) { $backbone2{$strain}++; }	# create a quick look up index for backbone2

open DATA, "<$datafile" or die "couldn't open $datafile : $!";
open DATA2, ">$hetdata" or die "couldn't open $hetdata : $!";

print DATA2 "strain\tfile\tmtl\tclade\tq40ps\tq40l\tq40het\tfps\tfl\tfhet\tfdps\tfdl\tfdhet\tlohl\n";


while (<DATA>) {
	if (/^(\S+)/) { if (exists $backbone2{$1}) { print DATA2 $_; } }
	if (/^NCYC/) { print DATA2 $_;  }
	if (/^SC5314/) { print DATA2 $_;  }

}
close DATA2;
close DATA;

open RCMD, ">Hetanalysis.Rcmd" or die "couldn't open Hetanalysis.Rcmd : $!";

print RCMD "

pdf(\"Rplots$datafile.pdf\")
rm(list=ls())
data<-read.table(\"$datafile\",header=T)
temp<-data[grep(\"CEC3678\",data\$strain,invert=T),]	# drop these 2 low quality strains
temp2<-temp[grep(\"CEC3638\",temp\$strain,invert=T),]	# drop these 2 low quality strains
data2<-temp2[grep(\"NCYC597\",temp2\$strain,invert=T),]	#drop NCYC597

attach(data2)
head(data2)
tail(data2)
length(strain)
sort(strain)

oakq40het<-q40het[grep(\"NCYC4\",strain)]
clinicq40het<-q40het[grep(\"NCYC4\",strain,invert=T)]		# this includes NCYC597 and SC5314

oakq40l<-q40l[grep(\"NCYC4\",strain)]
clinicq40l<-q40l[grep(\"NCYC4\",strain,invert=T)]
oakq40l
sort(clinicq40l)
summary(oakq40l)
summary(clinicq40l)
length(oakq40l)
length(clinicq40l)
wilcox.test(oakq40l,clinicq40l)				# not sig diff in length of high qual sequence between oak and clinic

ncycq40l<-q40l[grep(\"NCYC\",strain)]
otherq40l<-q40l[grep(\"NCYC\",strain,invert=T)]
wilcox.test(ncycq40l,otherq40l)
wilcox.test(oakq40l,otherq40l)
cor.test(fhet,q40l)
model<-lm(fdhet~q40l)
summary(model)


length(oakq40het)
length(clinicq40het)
summary(oakq40het)
summary(clinicq40het)
data2[grep(\"NCYC4\",strain),]
wilcox.test(oakq40het,clinicq40het)


oakfhet<-fhet[grep(\"NCYC4\",strain)]
clinicfhet<-fhet[grep(\"NCYC4\",strain,invert=T)]		# this includes NCYC597 and WT_SC5314


length(oakfhet)
length(clinicfhet)
summary(oakfhet)
summary(clinicfhet)
data2[grep(\"NCYC4\",strain),]
wilcox.test(oakfhet,clinicfhet)


oakfdhet<-fdhet[grep(\"NCYC4\",strain)]
clinicfdhet<-fdhet[grep(\"NCYC4\",strain,invert=T)]		# this includes NCYC597 and WT_SC5314

length(oakfdhet)
length(clinicfdhet)
summary(oakfdhet)
summary(clinicfdhet)
data2[grep(\"NCYC4\",strain),]
wilcox.test(oakfdhet,clinicfdhet)

o<-order(fdhet)
strain[o]
fdhet[o]


oakfdl<-fdl[grep(\"NCYC4\",strain)]
clinicfdl<-fdl[grep(\"NCYC4\",strain,invert=T)]		# this includes NCYC597 and WT_SC5314

length(oakfdl)
length(clinicfdl)
summary(oakfdl)
summary(clinicfdl)
data2[grep(\"NCYC4\",strain),]
wilcox.test(oakfdl,clinicfdl)

summary(fdhet\[clade==13\])

# Am I only seeing high heterozygosity for oaks because of the low heterozygosity of Clade 13?

no13data<-data2[clade!=13,]
clinicno13q40het<-no13data\$q40het[grep(\"NCYC4\",strain,invert=T)]		
wilcox.test(oakq40het,clinicno13q40het)

clinicno13fdhet<-no13data\$fdhet[grep(\"NCYC4\",strain,invert=T)]	
wilcox.test(oakfdhet,clinicno13fdhet)


# DO OAK STRAINS HAVE LESS LOH?	NO					

oaklohl<-lohl[grep(\"NCYC4\",strain)]
cliniclohl<-lohl[grep(\"NCYC4\",strain,invert=T)]		# this includes NCYC597 and WT_SC5314

length(oaklohl)
length(cliniclohl)
summary(oaklohl)
summary(cliniclohl)
data2[grep(\"NCYC4\",strain),]
wilcox.test(oaklohl,cliniclohl)


# ARE THERE DIFFERENCES BETWEEN CLADES? 
# YES, CLADE 13 is very different

nooakdata<-data2[grep(\"NCYC4\",strain,invert=T),]

# WITH LOH and REPEATS
model<-lm(nooakdata\$q40het~nooakdata\$clade)
summary(model)
model0<-lm(nooakdata\$q40het~1)
anova(model)
anova(model,model0)

model<-lm(nooakdata\$fdhet~nooakdata\$clade)
summary(model)
model0<-lm(nooakdata\$fdhet~1)
anova(model)
anova(model,model0)



# oak strains are more heterozygous than others in their clades

clade2<-clade
clade2[strain==\"NCYC4144\"]<-18
clade2[strain==\"NCYC4146\"]<-4
clade2[strain==\"NCYC4145\"]<-\"NC\"
model2<-lm(q40het~clade2)
summary(model2)

anova(model,model2)



# COMPARISON OF EACH OAK STRAIN WITH ITS OWN CLADE

q40het[strain==\"NCYC4146\"]
summary(q40het[clade==4&strain!=\"NCYC4146\"])
length(q40het[clade==4&strain!=\"NCYC4146\"])
sort(q40het[clade==4])
summary(fdl[clade==4&strain!=\"NCYC4146\"])
length(fdl[clade==4&strain!=\"NCYC4146\"])
summary(fdhet[clade==4&strain!=\"NCYC4146\"])
length(fdhet[clade==4&strain!=\"NCYC4146\"])



q40het[strain==\"NCYC4144\"]
summary(q40het[clade==18&strain!=\"NCYC4144\"])
length(q40het[clade==18&strain!=\"NCYC4144\"])
sort(q40het[clade==18])
summary(fdl[clade==18&strain!=\"NCYC4144\"])
length(fdl[clade==18&strain!=\"NCYC4144\"])
summary(fdhet[clade==18&strain!=\"NCYC4144\"])
length(fdhet[clade==18&strain!=\"NCYC4144\"])
fdhet[strain==\"NCYC4144\"]
sort(fdhet[clade==18])


q40het[strain==\"NCYC4145\"]
summary(q40het[clade==\"NC\"&strain!=\"NCYC4145\"])
length(q40het[clade==\"NC\"&strain!=\"NCYC4145\"])
sort(q40het[clade==\"NC\"])
summary(fdl[clade==\"NC\"&strain!=\"NCYC4145\"])
length(fdl[clade==\"NC\"&strain!=\"NCYC4145\"])
summary(fdhet[clade==\"NC\"&strain!=\"NCYC4145\"])
length(fdhet[clade==\"NC\"&strain!=\"NCYC4145\"])



# WITHOUT LOH and REPEATS
fmodel<-lm(fdhet~clade)
summary(fmodel)

# AFTER FILTERING oak strains are STILL more heterozygous than others in their clades

fmodel2<-lm(fdhet~clade2)
summary(fmodel2)

anova(fmodel,fmodel2)
anova(fmodel)

par(mfrow=c(2,2))
plot(fmodel)
plot(fmodel2)
#data2[14,]	# backbone analysis
#data2[31,]
#data2[15,]

data2[161,]	# all Ropars data analysis
data2[182,]
data2[98,]

hist(fdhet[clade==1],breaks=20)
hist(fdhet[clade==2],breaks=10)
hist(fdhet[clade==3],breaks=20)
hist(fdhet[clade==4],breaks=10)
hist(fdhet[clade==11],breaks=10)
hist(fdhet[clade==13],breaks=10)
hist(fdhet[clade==\"NC\"],breaks=10)


# FILTERED COMPARISON OF EACH CLADE WITH ITS OWN CLADE

fdhet[strain==\"NCYC4146\"]
summary(fdhet[clade==4])
length(fdhet[clade==4])
fdhet[clade==4]

fdhet[strain==\"NCYC4144\"]
summary(fdhet[clade==18])
length(fdhet[clade==18])
fdhet[clade==18]

fdhet[strain==\"NCYC4145\"]
summary(fdhet[clade==\"NC\"])
length(fdhet[clade==\"NC\"])
fdhet[clade==\"NC\"]



# EVEN after removing CLADE 13, clades are different (but less different than before) .. now oak strains do look the most different

attach(no13data)
model<-lm(q40het~clade)
summary(model)

tapply(q40het,clade,median)
tapply(q40het,clade,mean)

# EVEN after removing oak strains, clades are different 

noOakOr13data<-data2[data\$clade!=13&data\$clade!=\"oak\",]
attach(noOakOr13data)
model<-lm(q40het~clade)
summary(model)


# REMOVE CLADE 13 strains because they have lower variance than the rest

data2[data2\$clade==13,]
tapply(fdhet,clade,var)


data3<-data2[data2\$clade!=13,]
attach(data3)
length(data2\$strain)
length(data3\$strain)
length(fdhet)

fmodel<-lm(fdhet~clade)
summary(fmodel)

clade2<-clade
clade2[strain==\"NCYC4144\"]<-18
clade2[strain==\"NCYC4146\"]<-4
clade2[strain==\"NCYC4145\"]<-\"NC\"

fmodel2<-lm(fdhet~clade2)
summary(fmodel2)

anova(fmodel,fmodel2)

par(mfrow=c(2,2))
plot(fmodel)
plot(fmodel2)

# DO I GET THE SAME RESULT WITH A QUASIPOISSON GLM? this is not the right analysis because each heterozygous site is not independent

attach(data2)

fmodel<-glm(fdps~clade+offset(log(fdl)),quasipoisson)
summary(fmodel)

clade2<-clade
clade2[strain==\"NCYC4144\"]<-18
clade2[strain==\"NCYC4146\"]<-4
clade2[strain==\"NCYC4145\"]<-\"NC\"

fmodel2<-glm(fdps~clade2+offset(log(fdl)),quasipoisson)
summary(fmodel2)

anova(fmodel,fmodel2,test=\"Chi\")

# DO I GET THE SAME RESULT WITH A QUASIBINOMIAL GLM?

attach(data2)

fmodel<-glm(cbind(fdps,fdl)~clade,quasibinomial)
summary(fmodel)

clade2<-clade
clade2[strain==\"NCYC4144\"]<-18
clade2[strain==\"NCYC4146\"]<-4
clade2[strain==\"NCYC4145\"]<-\"NC\"

fmodel2<-glm(cbind(fdps,fdl)~clade2,quasibinomial)
summary(fmodel2)

anova(fmodel,fmodel2,test=\"Chi\")



";


close RCMD;

open SHOUT, ">runHetanalysis" or die "couldn't open runHetanalysis : $!";

printheader("runHetanalysis");

print SHOUT "`R CMD BATCH --no-save --no-restore Hetanalysis.Rcmd`\n";

close SHOUT;

print "\nNext step to run the R commands:\nqsub runHetanalysis\n";
print "Then results of analysis will be in in Hetanalysis.Rcmd.Rout, $hetdata and $table\n\n";


## REPLOT faChrompaint.pl plots with new colors and improve labels


open SHOUT, ">runRePaint" or die "couldn't open runRePaint : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runRePaint
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=64gb

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3

cd /lustre1/doudab/candidAnalysis/final2/paintWithOakGenomeM90

for strain in NCYC4144 NCYC4145 NCYC4146 NCYC597 SC5314 @backbone2
do
  faChrompaint.pl -p \$strain.nearest.tsv -c ../clades -C ../colors4 
  cat nearest\$strain.R | sed 's/mar=c(1,1,1,1)/mar=c(1,5.5,1,1),xpd=T/g' | sed \"s/main=\\\"\$strain.b/ylab='')\\ntext(-300000,50,labels=\\\"/g\" | sed \"s/.mfa\\\")/\\\",cex=3.5)/g\" > newnearest\$strain.R
 


  R CMD BATCH --no-restore newnearest\$strain.R
  mv \$strain.nearest.pdf \${strain}"."newnearest4.pdf
done


";

close SHOUT;

print "\nNext step:\nqsub runRePaint\n\n";

################
## run unique.pl on sapelo2
###############

open SHOUT, ">runUnique3" or die "couldn't open runUnique3 : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q bensasson_q
#PBS -N runUnique3
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -l mem=128gb

#PBS -M doudab\@uga.edu 
#PBS -m ae


cd /lustre1/doudab/candidAnalysis/final2/

for chr in @chr
do
 
 fastaLC2n.pl -i b\$chr.mfa -A -o amb2n\$chr.mfa
 perl al2variable.pl -a amb2n\$chr.mfa -g -o var\$chr
 perl unique.pl -i  var\$chr.gfa -s -o  uniqvar\$chr -R
done

";

close SHOUT;

print "\nNext step:\nqsub runUnique3\n\n";




#################
## SUBROUTINES ##
#################


sub printphyloheader {

	my $runname = shift;

	print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N $runname
#PBS -l nodes=1:ppn=16:dbnode
#PBS -l walltime=168:00:00
#PBS -l mem=16gb

#PBS -M doudab\@uga.edu 
#PBS -m ae


cd $dir

module load raxml/8.1.20
module load R/3.2.3
";

	
}

sub printheader {

	my $runname = shift;

	print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N $runname
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=8gb

#PBS -M doudab\@uga.edu 
#PBS -m ae


cd $dir

module load R/3.2.3
";

}


# SUBROUTINE for drawing trees in R

sub printRcmds4fig {

	my $chr = shift;

	my %nodes = (
		"chr1" => 35,
		"chr2" => 34,
		"chr3" => 43,
		"chr4" => 34,
		"chr5" => 33,
		"chr6" => 39,
		"chr7" => 34,
		"chrR" => 45,
#		"chr1chr2chr3chr4chr5chr6chr7chrR" => 36,
		"chr1chr2chr3chr4chr5chr6chr7chrR" => 31,
		"MLSThir" => 35,
		"chr5_150000_250000" => 39,
		"chrR_2000000_2100000" => 46,
		"hirncycanimaldsts" => 70,
		"genome" => 50,
		"genome2" => 98,

	);

	if ($chr eq "hirncycanimaldsts.mfa.cons") { print RCMD "t<-read.tree('$chr.tree')\n"; }
	else {
		print RCMD "

t<-read.tree('raxml$chr/RAxML_bipartitions.$chr')

#plot(t)					# Figure out node names for re-rooting
#nodelabels(as.character((1:length(t\$node))+length(t\$tip)))
";
	}

	if (defined $nodes{$chr}) { print RCMD "rt<-root(t,node=$nodes{$chr})\nt<-rt\n"; } 	# re-root the tree

	print RCMD "

clade1col<-\"".$colors{"1"}."\"
clade2col<-\"".$colors{"2"}."\"			
clade3col<-\"".$colors{"3"}."\"
clade4col<-\"".$colors{"4"}."\" 			# light blue
clade8col<-\"".$colors{"8"}."\"		
clade9col<-\"".$colors{"9"}."\"
clade10col<-\"".$colors{"10"}."\"
clade11col<-\"".$colors{"11"}."\"	
clade12col<-\"".$colors{"12"}."\"
clade13col<-\"".$colors{"13"}."\"
clade16col<-\"".$colors{"16"}."\"
clade18col<-\"".$colors{"18"}."\"
cladeAcol<-\"".$colors{"A"}."\"
cladeBcol<-\"".$colors{"B"}."\"
cladeCcol<-\"".$colors{"C"}."\"
cladeDcol<-\"".$colors{"D"}."\"
cladeEcol<-\"".$colors{"E"}."\"

unknowncol<-\"".$colors{"NC"}."\"

wildcol<-\"".$colors{"oak"}."\"			# green
#ncyccol<-\"".$colors{""}."\"
tc<-rep(\"$ec\",length(t\$tip)) 	# default tip color = unknown
ew<-rep($ew,length(t\$edge[,2]))		# default edge width = 1
ec<-rep(\"$ec\",length(t\$edge[,2]))	# default edge color = black

pew<-$pew					# population edge width


tc[grep(\"NCYC4\",t\$tip)]<-wildcol	


tc[grep(\"CEC3682\",t\$tip)]<-clade1col
tc[grep(\"CEC3623\",t\$tip)]<-clade1col
tc[grep(\"CEC3560\",t\$tip)]<-clade1col
tc[grep(\"SC5314\",t\$tip)]<-clade1col

tc[grep(\"CEC4946\",t\$tip)]<-clade2col
tc[grep(\"CEC3549\",t\$tip)]<-clade2col
tc[grep(\"CEC3531\",t\$tip)]<-clade2col

tc[grep(\"CEC3637\",t\$tip)]<-clade3col
tc[grep(\"CEC3681\",t\$tip)]<-clade3col
tc[grep(\"CEC3557\",t\$tip)]<-clade3col

tc[grep(\"CEC3676\",t\$tip)]<-clade4col
tc[grep(\"CEC3607\",t\$tip)]<-clade4col
tc[grep(\"CEC3674\",t\$tip)]<-clade4col


tc[grep(\"CEC2023\",t\$tip)]<-clade8col
tc[grep(\"CEC3622\",t\$tip)]<-clade8col
tc[grep(\"CEC3634\",t\$tip)]<-clade8col

tc[grep(\"CEC2018\",t\$tip)]<-clade9col
tc[grep(\"CEC3533\",t\$tip)]<-clade9col
tc[grep(\"CEC5026\",t\$tip)]<-clade9col

tc[grep(\"CEC3616\",t\$tip)]<-clade10col
tc[grep(\"CEC3706\",t\$tip)]<-clade10col
tc[grep(\"CEC3711\",t\$tip)]<-clade10col

tc[grep(\"CEC4495\",t\$tip)]<-clade11col
tc[grep(\"CEC3618\",t\$tip)]<-clade11col
tc[grep(\"CEC3704\",t\$tip)]<-clade11col

tc[grep(\"CEC3680\",t\$tip)]<-clade12col
tc[grep(\"CEC2019\",t\$tip)]<-clade12col
tc[grep(\"CEC3537\",t\$tip)]<-clade12col

tc[grep(\"CEC4877\",t\$tip)]<-clade13col
tc[grep(\"CEC5020\",t\$tip)]<-clade13col
tc[grep(\"CEC5019\",t\$tip)]<-clade13col

tc[grep(\"CEC3550\",t\$tip)]<-clade16col
tc[grep(\"CEC3600\",t\$tip)]<-clade16col
tc[grep(\"CEC3664\",t\$tip)]<-clade16col

tc[grep(\"CEC3786\",t\$tip)]<-clade18col
tc[grep(\"CEC2871\",t\$tip)]<-clade18col
tc[grep(\"CEC2876\",t\$tip)]<-clade18col

tc[grep(\"CEC3548\",t\$tip)]<-cladeAcol
tc[grep(\"CEC3715\",t\$tip)]<-cladeAcol
tc[grep(\"CEC4038\",t\$tip)]<-cladeAcol

tc[grep(\"CEC3708\",t\$tip)]<-cladeBcol
tc[grep(\"CEC4254\",t\$tip)]<-cladeBcol

tc[grep(\"CEC3661\",t\$tip)]<-cladeCcol
tc[grep(\"CEC4039\",t\$tip)]<-cladeCcol

tc[grep(\"CEC3707\",t\$tip)]<-cladeDcol
tc[grep(\"CEC4498\",t\$tip)]<-cladeDcol

tc[grep(\"CEC3712\",t\$tip)]<-cladeEcol
tc[grep(\"CEC4486\",t\$tip)]<-cladeEcol

tc[grep(\"CEC3555\",t\$tip)]<-unknowncol
tc[grep(\"CEC3613\",t\$tip)]<-unknowncol
tc[grep(\"CEC3619\",t\$tip)]<-unknowncol
tc[grep(\"CEC3663\",t\$tip)]<-unknowncol
tc[grep(\"CEC3665\",t\$tip)]<-unknowncol
tc[grep(\"CEC3671\",t\$tip)]<-unknowncol
tc[grep(\"CEC4261\",t\$tip)]<-unknowncol
tc[grep(\"CEC4481\",t\$tip)]<-unknowncol
tc[grep(\"CEC4485\",t\$tip)]<-unknowncol
tc[grep(\"CEC723\",t\$tip)]<-unknowncol

    
							# color the clades by clade
							# function from http://blog.phytools.org/2012/01/function-to-get-descendant-node-numbers.html
getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree\$edge[which(tree\$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree\$tip))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}
";

	
	foreach my $clade (sort keys %node) {
		if ($clade =~ /[BCDE][12]/) { next; }
		print RCMD "
cladenode$clade<-$node{$clade}
d<-getDescendants(t,node=cladenode$clade)
for (x in d) { ec[t\$edge[,2]==x]<-clade$clade"."col; ew[t\$edge[,2]==x]<-pew }
ec[t\$edge[,2]==cladenode$clade]<-clade$clade"."col		# color the edge leading to this clade
ew[t\$edge[,2]==cladenode$clade]<-pew
	
";
	}

	print RCMD "

ec[t\$edge[,2]==".$node{"B1"}."]<-cladeBcol		# color in edges for clades with only 2 strains
ec[t\$edge[,2]==".$node{"B2"}."]<-cladeBcol
ec[t\$edge[,2]==".$node{"C1"}."]<-cladeCcol
ec[t\$edge[,2]==".$node{"C2"}."]<-cladeCcol
ec[t\$edge[,2]==".$node{"D1"}."]<-cladeDcol
ec[t\$edge[,2]==".$node{"D2"}."]<-cladeDcol
ec[t\$edge[,2]==".$node{"E1"}."]<-cladeEcol
ec[t\$edge[,2]==".$node{"E2"}."]<-cladeEcol



plot(t,edge.color=ec,tip.color=tc,edge.width=ew,cex=0.8)		# colour tips by clade

bs<-as.numeric(t\$node)				# show only bootstraps >= 70%
bs[bs<70]<-NA
#nodelabels(bs,frame=\"none\",cex=0.5,bg=\"white\",adj = c(1.1,-0.2))
add.scale.bar(0.002,1)

#nodelabels(bs,frame=\"none\",bg=\"white\",adj = c(1.5,-0.3))
nodelabels(bs,frame=\"none\",bg=\"white\",adj = c(1.1,-0.2),cex=0.8,col=\"$ec\")

";

}



