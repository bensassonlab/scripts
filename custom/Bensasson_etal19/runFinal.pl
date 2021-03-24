#!/usr/bin/perl

use warnings;
use strict;

my $program = 'runFinal.pl';					#name of script

my @NCYC = qw(NCYC597 NCYC4144 NCYC4145 NCYC4146);

my $ncycdir = "/lustre1/doudab/candidAnalysis/NCYC/";
my $skip = "yes";
my $refseqs = "C_albicans_SC5314_version_A22-s07-m01-r18_chromosomes.fasta";
my $refdir = "refAandB";
my $hethist = "hethist";
my $maxLOH = 0.001;
my $withinCladeMax = 0.0006604028;
my $bootstraps = 1000;
my $tiff = "yes";
my $datafile = "data.LOHnoR.tsv"; # OUTPUT TABLE OF DATA FOR EACH STRAIN FROM HETEROZYGOSITY ANALYSIS
my $table = "Table2.tsv";
my $chrlengths = "c(3187324+1000000,2231865,1799223,1602610,1190833,1032969,949245,2285932)"; # add some length to chr1 for the margins

my $tableDH = "Table2ii.tsv";
my $stable = "TableS1.tsv";

# NOTE: faChrompaint.pl requires 2 files: 'colors' and 'clades'


# POSITIONS OF MLST LOCI ON REFERENCE GENOME

my $mlstpos = "
Locus	Allele	Length	Contig	Start position	End position
AAT1a	8	373	chr2	1078561	1078933
ACC1	3	407	chrR	158446	158852
ADP1	2	443	chrR	1247585	1248027
MPIb	9	375	chr2	1972847	1973221
SYA1	2	391	chr6	806597	806987
VPS13	24	403	chr4	1346953	1347355
ZWF1b	12	491	chr1	1953391	1953881
";
my %loci = (
	AAT1a => "revcomp",
	ACC1 => "fwd",
	ADP1 => "revcomp",
	MPIb => "fwd",
	SYA1 => "fwd",
	VPS13 => "fwd",
	ZWF1b => "revcomp",
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

my %clades = (
	"NCYC4144" => "?",
	"NCYC4145" => "?",
	"NCYC4146" => 4,
	"WT_SC5314" => 1,
	"12C" => 1,
	"L26" => 1,
	"P37005" => 1,
	"19F" => 1,
	"P78048" => 1,
	"P37037" => 1,
	"P37039" => 1,
	"P57072" => 2,
	"P78042" => 3,
	"P34048" => 3,
	"P87" => 4,
	"GC75" => 4,
	"P75016" => 4,
	"P75063" => 4,
	"P94015" => 6,
	"P60002" => 8,
	"P75010" => 11,
	"P57055" => 3,
	"P76055" => 2,
	"P76067" => 2,
	"SC5314" => 1,
	"1AA" => 1,
);




# GET START FILES

unless ($skip eq "yes") {

`cp ../P78042/SC5314_A22edited.mfa .`;
`mkdir 12C_PRJNA75209`;
`mkdir 19F_PRJNA75221`;
`mkdir GC75_PRJNA75223`;
`mkdir L26_PRJNA75211`;
`mkdir P34048_PRJNA75229`;
`mkdir P37005_PRJNA75217`;
`mkdir P37037_PRJNA75231`;
`mkdir P37039_PRJNA75233`;
`mkdir P57055_PRJNA75239`;
`mkdir P57072_PRJNA75227`;
`mkdir P60002_PRJNA75219`;
`mkdir P75010_PRJNA75235`;
`mkdir P75016_PRJNA75237`;
`mkdir P75063_PRJNA75241`;
`mkdir P76055_PRJNA75243`;
`mkdir P76067_PRJNA75245`;
`mkdir P78048_PRJNA75225`;
`mkdir P87_PRJNA75215`;
`mkdir P94015_PRJNA75213`;
`cp ../Hirakawa_etal14/SC5314annotations.gff .`;
`mkdir SC5314`;
`mkdir 1AA`;
`ln -s /lustre1/doudab/candidAnalysis/SC5314/mapSRR850115toref/*_1.fastq.gz 1AA/`;
`ln -s /lustre1/doudab/candidAnalysis/SC5314/mapSRR850115toref/*_2.fastq.gz 1AA/`;
`ln -s /lustre1/doudab/candidAnalysis/SC5314/mapSRR850113toref/*_1.fastq.gz SC5314/`;
`ln -s /lustre1/doudab/candidAnalysis/SC5314/mapSRR850113toref/*_2.fastq.gz SC5314/`;
`mkdir P78042_PRJNA75247`;
`ln -s /lustre1/doudab/candidAnalysis/P78042/P78042_PRJNA75247/*_1.fastq.gz P78042_PRJNA75247/`;
`ln -s /lustre1/doudab/candidAnalysis/P78042/P78042_PRJNA75247/*_2.fastq.gz P78042_PRJNA75247/`;
}

# SET UP SCRIPTS FOR TRIMMING READS

foreach my $strain (@NCYC) {

	my $afile = $strain."adapters.fa";

# GET MORE START FILES FOR EACH STRAIN

	unless ($skip eq "yes") {
		`mkdir $strain`;
		`ln -s $ncycdir/$strain/*_R1.fastq.gz $strain/.`;
		`ln -s $ncycdir/$strain/*_R2.fastq.gz $strain/.`;
		`ln -s $ncycdir/$strain/trimmed/$afile $strain/.`;
	}

	
	open OUT, ">runtrim$strain.sh" or die "couldn't open runtrim$strain.sh : $!";

	print OUT "
#PBS -S /bin/bash
#PBS -q batch
#PBS -N trim$strain
#PBS -l nodes=1:ppn=4:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=64gb		# really need this much

#PBS -M doudab\@uga.edu 
#PBS -m ae

cd /lustre1/doudab/candidAnalysis/final/$strain

module load trimmomatic/0.33

for sample in *_R1.fastq.gz
do
 i=`basename \$sample _R1.fastq.gz`
 input1=\"\${i}_R1.fastq.gz\"
 input2=\"\${i}_R2.fastq.gz\"
 output1p=\"\${i}_R1_trimmed_paired.fastq.gz\"
 output1u=\"\${i}_R1_trimmed_unpaired.fastq.gz\"
 output2p=\"\${i}_R2_trimmed_paired.fastq.gz\"
 output2u=\"\${i}_R2_trimmed_unpaired.fastq.gz\"
 

 time java -jar /usr/local/apps/trimmomatic/0.33/trimmomatic-0.33.jar PE -threads 4 \$input1 \$input2 \$output1p \$output1u \$output2p \$output2u ILLUMINACLIP:$afile:2:30:10 

done

end
";

	close OUT;
}

print "\nNext steps for trimming reads:\n";

foreach my $strain (@NCYC) {
	print "qsub runtrim$strain.sh\n";
}

print "\n";

# SET UP SCRIPTS FOR MAPPING AND BASECALLING

opendir (DIR, ".") or die "couldn't open . : $!";

my (%strains,@dirs);	
while (defined (my $name = readdir(DIR))) {

	if ($name =~ /^([a-z0-9]+)_([a-z0-9]+)$/i) { 
		print "$name\tstrain:$1 accession:$2\n"; 	
									# Skip strains that already have a binary bam file in their directory
		if (-B "$name/$1.bam") { print "Skipping because already mapped: $name/$1.bam\n"; next; }
		push (@dirs, $name);
		$strains{$name} = $1;
		unless ($skip eq "yes") { `ln -s /lustre1/doudab/candidAnalysis/Hirakawa_etal14/$name/SRR*fastq.gz $name/.`; }
	}
	elsif ($name =~ /^(NCYC\d+)$/) {
		print "$name\tstrain:$1\n"; 	
									# Skip strains that already have a binary bam file in their directory
		if (-B "$name/$1.bam") { print "Skipping because already mapped: $name/$1.bam\n"; next; }
		push (@dirs, $name);
		$strains{$name} = $1;
		
	}
	elsif (($name eq "SC5314") || ($name eq "1AA")) {
		print "$name\tstrain:$name\n"; 	
									# Skip strains that already have a binary bam file in their directory
		if (-B "$name/$name.bam") { print "Skipping because already mapped: $name/$name.bam\n"; next; }
		push (@dirs, $name);
		$strains{$name} = $name;

	}
}		

print "Next steps for mapping reads, basecalling, fasta consensus and aligning:\n";

foreach my $dir (@dirs) {

	open SHOUT, ">runMapCall$strains{$dir}" or die "couldn't open runMapCall$strains{$dir} : $!";

	print SHOUT "
#PBS -S /bin/bash
#PBS -q batch
#PBS -N map$strains{$dir}
#PBS -l nodes=1:ppn=1:dbnode		# reduce nodes because samtools takes ages and only uses 1 processor
#PBS -l walltime=168:00:00
#PBS -l mem=16gb			# job history on these does not suggest very high memory

#PBS -M doudab\@uga.edu 
#PBS -m ae


cd /lustre1/doudab/candidAnalysis/final/$dir

module load bwa/0.7.10
module load samtools/1.2
module load freebayes/1.0.1
module load R/3.2.3
module load seqtk/1.0-r82-dirty

## MAP THE HIRAKAWA AND SC5314 DATA

for sample in *_1.fastq.gz
do
 i=`basename \$sample _1.fastq.gz`
 input1=\"\${i}_1.fastq.gz\"
 input2=\"\${i}_2.fastq.gz\"
 output=\"\${i}\"
 
 bwa mem -t 1 /lustre1/doudab/candidAnalysis/final/SC5314_A22edited.mfa \$input1 \$input2 | samtools view -b - | samtools sort - \$output

done

## MAP THE NCYC DATA

for sample in *_R1_trimmed_paired.fastq.gz
do
 i=`basename \$sample _R1_trimmed_paired.fastq.gz`
 input1=\"\${i}_R1_trimmed_paired.fastq.gz\"
 input2=\"\${i}_R2_trimmed_paired.fastq.gz\"
 output=\"\${i}\" 

 bwa mem -t 1 /lustre1/doudab/candidAnalysis/final/SC5314_A22edited.mfa \$input1 \$input2 | samtools view -b - | samtools sort - \$output

done

";

	if ($dir =~ /NCYC/) {
		print SHOUT "
#rm $strains{$dir}.bam
samtools merge $strains{$dir}.bam Run*.bam
rm Run*.bam
";
	}
	elsif ($dir eq "1AA")  { print SHOUT "mv SRR*.bam 1AA.bam\n"; }		# only 1 pair of reads for each of these, therefore no need for merging bam files
	elsif ($dir eq "SC5314") { print SHOUT "mv SRR*.bam SC5314.bam\n"; }
	else {

		print SHOUT "
#rm $strains{$dir}.bam
samtools merge $strains{$dir}.bam SRR*.bam
rm SRR*.bam
";
	}
	print SHOUT "

samtools mpileup -I -d 10000 -uf /lustre1/doudab/candidAnalysis/final/SC5314_A22edited.mfa $strains{$dir}.bam | bcftools call -c > $strains{$dir}.I.vcf

# create a link for subsequent vcf analysis in the main directory
ln -s $strains{$dir}.I.vcf ..

# VISUALIZE THE VCF FILE
ln -s ../SC5314annotations.gff .
vcf2allelePlot.pl -i $strains{$dir}.I.vcf -g SC5314annotations.gff -h -m

# CONVERT VCF TO FASTA
vcfutils.pl vcf2fq $strains{$dir}.I.vcf > $strains{$dir}.I.fq
 seqtk seq -q 40 -A $strains{$dir}.I.fq >$strains{$dir}.I.fa

# CREATE ALIGNMENTS FOR EACH CHROMOSOME
for chr in chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chrR
do
  faChoose.pl -i $strains{$dir}.I.fa -n \$chr -o $strains{$dir}\$chr.I.fa 
  cat $strains{$dir}\$chr.I.fa | sed s/\$chr/$strains{$dir}/ > temp
  mv temp $strains{$dir}\$chr.I.fa
done
echo \"created 8 $strains{$dir} files\"
";

	print "qsub runMapCall$strains{$dir}\n";

	close SHOUT;	

}


## CREATE LOH FILTER AND MORE VCF2ALLELES PLOTS

opendir (DIR, ".") or die "couldn't open . : $!";

my @vcffiles;						# collect the names of data directories, strains and their study accessions
while (defined (my $name = readdir(DIR))) {

	if ($name =~ /^(\S+)\.I\.vcf$/i) { 
		print "$name\tstrain:$1\n"; 	
		$strains{$name} = $1;
		push (@vcffiles, $name);
	}
}		





# estimate heterozygosity for every strain after filtering out LOH regions
print "# The following scripts will visualize the vcf plots and\n";
print "# estimate heterozygosity for every strain after filtering out LOH regions\n";

print "\n# Next step:\n\n";

foreach my $file (@vcffiles) { 



	open SHOUT, ">runVcfPlots$strains{$file}" or die "couldn't open runVcfPlots$strains{$file} : $!";
	print SHOUT "
#PBS -S /bin/bash
#PBS -q batch
#PBS -N vcf2Plots$strains{$file}
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=16gb		# these jobs usually only require 5GB of memory, and increasing memory increases wait in the queue

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3		# vcf2allelePlot.pl uses R



cd /lustre1/doudab/candidAnalysis/final

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

#########################
### DRAW ALLELE PLOT FIGURES FOR MANUSCRIPT
##########################

# create allele plots with LOH for NCYC strains for Fig 1 and for other strains

my @chr = qw(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chrR);


open SHOUT, ">runallelePlotFigs" or die "couldn't open runallelePlotFigs : $!";
print SHOUT "
#PBS -S /bin/bash
#PBS -q batch
#PBS -N allelePlotFigs
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=16gb		# these jobs usually only require 5GB of memory, and increasing memory increases wait in the queue

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3		# vcf2allelePlot.pl uses R



cd /lustre1/doudab/candidAnalysis/final
";

my %ncyc = (
	"NCYC4144" => 1,
	"NCYC4145" => 1,
	"NCYC4146" => 1,
#	"NCYC597" =>1,
);
my %nonhirakawa = (
	"SC5314" => 1,
	"1AA" => 1,
);

opendir (DIR, ".") or die "couldn't open . : $!";

my (@strains,@callfiles);					# collect the names of data directories, strains and their study accessions
while (defined (my $name = readdir(DIR))) {

	if ($name =~ /^(\S+)\.I\.calls$/i) { 
		print "$name\tstrain:$1\n"; 	
		$strains{$name} = $1;
		push (@callfiles, $name);
		push (@strains,$1);
	}

	if ($name =~ /^(\S+)\.LOH\.I\.gff$/i) { 
		print "$name\tstrain:$1\n"; 	
		$strains{$name} = $1;
	}

}		

open RCMD, ">gwdCallPlot.Rcmd" or die "couldn't open gwdCallPlot.Rcmd : $!";

print RCMD "
#pdf(\"NCYCplots.pdf\")
#tiff(\"NCYCplots.tiff\")	# does not work
#jpeg(\"NCYCplots.jpeg\",width=2000,height=2000)

#bitmap(\"NCYCplots.tiff\", height = 19, width = 19, units = 'cm', type = \"tiff24nc\", res = 600)
#bitmap(\"NCYCplots.tiff\", height = 23, width = 19, units = 'cm', type = \"tiff24nc\", res = 600)
bitmap(\"NCYCplots.tiff\", height = 12, width = 19, units = 'cm', type = \"tiff24nc\", res = 600)


#par(mfrow=c(4,8))
#layout(rbind(1:8,9:16,17:24,25:32),widths=$chrlengths) 			# 4 strains per page
#layout(rbind(1:8,9:16,17:24,25:32,33:40,41:48),widths=$chrlengths)		# 6 strains per page
layout(rbind(1:8,9:16,17:24),widths=$chrlengths)			å	# 3 strains per page



options(scipen=999)
";

my $seenNCYC;
my (%atype,%astart,%aend,%seentype);		# chr will be the keys in these hashes of arrays (these are needed in the subroutines below and need to be global

foreach my $file (@callfiles) {						# read annotations for all files

	my $gffile = "$strains{$file}.LOH.I.gff";

	if (defined $gffile) { 			# Read annotations if they exist
		open GFF, "<$gffile" or die "couldnt open $gffile : $!"; 
	#	print "Reading in LOH summary from $gffile\n";
		while (<GFF>) {
			if (/^(\S+)\s+\S+\s+(LOH)\s+(\d+)\s+(\d+)/m) {
				$seentype{$2}++;
				push (@{$atype{"$file$1"}},$2);
				push (@{$astart{"$file$1"}},$3);
				push (@{$aend{"$file$1"}},$4);			
			}				
		}
		unless (defined $seentype{"LOH"}) { warn "ERROR: $gffile is empty or has an unrecognized format!\n"; }

	}
}

foreach my $file (sort @callfiles) {
	if (defined $ncyc{$strains{$file}}) { printstrain($file); }	# Print NCYC strains
}

foreach my $file (@callfiles) {						# Print SC5314-derived strains
#	if (defined $nonhirakawa{$strains{$file}}) { if ($file =~ /SC5314/) { printstrain($file); }} # comment this out to put SC5314-derived strains in the supplement
}

foreach my $file (@callfiles) {						# Print SC5314-derived strains
#	if (defined $nonhirakawa{$strains{$file}}) { if ($file =~ /1AA/) { printstrain($file); }} # comment this out to put SC5314-derived strains in the supplement
}



foreach my $file (sort @callfiles) {						# Print Hirakawa strains
#	unless ((defined $ncyc{$strains{$file}}) || (defined $nonhirakawa{$strains{$file}})) {  
	unless ((defined $ncyc{$strains{$file}})) {  	
		unless ($seenNCYC++) { 
			print RCMD "
#pdf(\"nonNCYCplots.pdf\")

#bitmap(\"nonNCYCplots.tiff\", height = 19, width = 19, units = 'cm', type = \"tiff24nc\", res = 300)
bitmap(\"nonNCYCplots.tiff\", height = 23, width = 19, units = 'cm', type = \"tiff24nc\", res = 600)


#par(mfrow=c(4,8),cex.lab=0.7,cex.axis=0.7,cex.main=0.7)
#layout(rbind(1:8,9:16,17:24,25:32),widths=$chrlengths)			# 4 strains per page
layout(rbind(1:8,9:16,17:24,25:32,33:40,41:48),widths=$chrlengths)	# 6 strains per page

options(scipen=999)
"; 
		}
		printstrain($file); 
	}
}

foreach my $file (@callfiles) {						# Print SC5314-derived strains
#	if (defined $nonhirakawa{$strains{$file}}) { printstrain($file); } # uncomment this to put SC5314-derived strains in the supplement
}


close RCMD;

print SHOUT "`R CMD BATCH --no-save --no-restore gwdCallPlot.Rcmd`\n";
print SHOUT "`tiff2pdf NCYCplots.tiff -o NCYCplots.pdf`\n";
print SHOUT "`tiff2pdf nonNCYCplots.tiff -o nonNCYCplots.pdf`\n";
close SHOUT;

print "\nNext step:\nqsub runallelePlotFigs\n\n";



# CREATE CHROMOSOME AND GENOMEWIDE ALIGNMENTS


open SHOUT, ">runPhylo" or die "couldn't open runPhylo : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runPhylo
#PBS -l nodes=1:ppn=4:dbnode
#PBS -l walltime=168:00:00
#PBS -l mem=64gb

#PBS -M doudab\@uga.edu 
#PBS -m ae


cd  /lustre1/doudab/candidAnalysis/final/

module load raxml/8.1.20
module load R/3.2.3
";

my ($allchrfiles,$allname);

open RCMD, ">showTree.R" or die "couldn't open showTree.R : $!";


if ($tiff eq "yes") { print RCMD "bitmap(\"raxmltrees.tiff\", height = 7, width = 7, units = 'in', type = \"tiff24nc\", res = 300)\n"; }
else { print RCMD "pdf('raxmltrees.pdf')\n"; }

print RCMD "

library(ape)
par(mfrow=c(2,2),mar=c(1,1,1,1),cex=0.7,xpd=T)
";

my $run = "short";
#my $run = "long";	# uncheck for a full/long run



if ($run eq "long") {
# prep haplotypes A and B for SC5314
#
	foreach my $chr (@chr) {

		print SHOUT "
faChoose.pl -i */$refseqs -n $chr"."A -o $refdir/SC5314_A$chr.I.fa 
cat $refdir/SC5314_A$chr.I.fa | sed 's/.*chr.*/>SC5314_A/g' > temp
mv temp $refdir/SC5314_A$chr.I.fa


";
	}
}

# Prep whole chromosome alignment files and run raxml on each one
my $seenlegend;
foreach my $chr (@chr) {

	if ($run eq "long") {


		$allchrfiles .= "$chr.mfa "; $allname .= $chr;

		`mkdir raxml$chr/`;

		print SHOUT "
cat */*$chr.I.fa > $chr.mfa
fastaLC2n.pl -i $chr.mfa -o temp
alcat.pl -i temp -o $chr.mfa
fa2phylip.pl -i $chr.mfa > $chr.phy

rm raxml$chr/*chr*
raxmlHPC-PTHREADS-AVX -T 4 -f a -x 12345 -p 12345 -m GTRGAMMA -N $bootstraps -s $chr.phy -n $chr -w /lustre1/doudab/candidAnalysis/final/raxml$chr/
";
	}
	unless ($chr =~/chr/) { die " unrecognized chromosome name: [$chr]\n"; }
	printRcmds4fig($chr);


	print RCMD "legend(\"bottomleft\", \"$chr\", bty=\"n\",cex=1.5,inset=0.05)\n";
	unless ($seenlegend++) { 
		print RCMD "legend(\"bottomright\",c(\"oak\",\"clade1\",\"clade2\",\"clade3\",\"clade4\",\"clade6\",\"clade8\",\"clade11\",\"unknown\"),text.col=c(wildcol,clade1col,clade2col,clade3col,clade4col,clade6col,clade8col,clade11col,ncyccol),bty=\"n\",cex=1.5)\n"; 
	}

}

## PRINT A LARGE LEGEND FOR FIG S7

print RCMD "
par(mfrow=c(1,1),mar=c(1,1,1,1))
plot(1,1,type='n',yaxt='n',xlab='',xaxt='n',bty='n')
legend(\"center\",c(\"oak\",\"clade1\",\"clade2\",\"clade3\",\"clade4\",\"clade6\",\"clade8\",\"clade11\",\"unknown\"),text.col=c(wildcol,clade1col,clade2col,clade3col,clade4col,clade6col,clade8col,clade11col,ncyccol),bty=\"n\",cex=4)
";


close RCMD;


print SHOUT "
`R CMD BATCH --no-save --no-restore showTree.R`
";

if ($tiff eq "yes") { print SHOUT "`tiff2pdf raxmltrees.tiff -o raxmltrees.pdf`\n"; }

close SHOUT;

print "\nNext step:\nqsub runPhylo\n\n";


# ANALYSIS OF MLST REGIONS

# READ IN POSITIONS OF EACH MLST ON NCYC AND HIRAKAWA GENOMES FROM REFERENCE (ncbiSC5314_A22.tsv)

$allname = "MLSThir";
my $dir = "/lustre1/doudab/candidAnalysis/final";

`mkdir raxmlMLSThir`;

open SHOUT, ">run$allname" or die "couldn't open run$allname : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N run$allname
#PBS -l nodes=1:ppn=4:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=64gb

#PBS -M doudab\@uga.edu 
#PBS -m ae


cd /lustre1/doudab/candidAnalysis/final/

module load raxml/8.1.20
module load R/3.2.3
";

my (@loci,%chr,%starts,%ends);

while ($mlstpos =~ /(\S+)\s+\d+\s+\d+\s+(\S+)\s+(\d+)\s+(\d+)/gs) { 
		push(@loci, $1);
		$chr{$1} = $2;
		$starts{$1} = $3;
		$ends{$1} = $4;
		print "MLST locus: $1 chr: $chr{$1} $starts{$1}..$ends{$1}\n";
}


foreach my $strain (@strains) {

	foreach my $locus (@loci) {
		if ($loci{$locus} eq "revcomp") {
			print SHOUT "faChooseSubseq.pl -i */$strain.I.fa -n $chr{$locus} -s '$starts{$locus}..$ends{$locus}' -o $locus.$strain.fa -r\n";
		}
		else {
			print SHOUT "faChooseSubseq.pl -i */$strain.I.fa -n $chr{$locus} -s '$starts{$locus}..$ends{$locus}' -o $locus.$strain.fa\n";
		}
#		print "\t$locus.$file\n";
		print SHOUT "cat $locus.$strain.fa | sed 's/chr.*/$strain/' > temp\n";	# rename sequence so only strain name is given (needed this way for alcat.pl later).
		print SHOUT "mv temp $locus.$strain.fa\n";
	}


}

# Cat-align into multiple alignments for each chromosome

print "\n";
my @mfafiles;
foreach my $locus (@loci) { 
	print SHOUT "cat $locus.*.fa > $locus.mfa\n"; 
	print SHOUT "fastaLC2n.pl -i $locus.mfa -o $locus.LC2n.mfa\n";
	print SHOUT "fa2phylip.pl -i $locus.LC2n.mfa > $locus.LC2n.phy\n";
	print SHOUT "rm $locus.*.fa\n";							# tidy up
	print SHOUT "rm $locus.mfa\n";							# tidy up
	print "$locus alignment is in $locus.LC2n.mfa and $locus.LC2n.phy\n";
	push(@mfafiles,"$locus.LC2n.mfa");
}

# Concatenate into a single alignment

print SHOUT "alcat.pl -i '@mfafiles' -o $allname.mfa\n";
print SHOUT "fa2phylip.pl -i $allname.mfa > $allname.phy\n";

print SHOUT "

rm raxmlMLSThir/RAxML*
raxmlHPC-PTHREADS-AVX -T 4 -f a -x 12345 -p 12345 -m GTRGAMMA -N $bootstraps -s $allname.phy -n $allname -w $dir/raxmlMLSThir/

";

open RCMD, ">raxmlMLSThir/showTree$allname.R" or die "couldn't open raxmlMLSThir/showTree$allname.R : $!";

print RCMD "

library(ape)
pdf('raxml$allname.pdf')

";


printRcmds4fig($allname);

print RCMD "
legend(\"bottomright\",c(\"oak\",\"clade1\",\"clade2\",\"clade3\",\"clade4\",\"clade6\",\"clade8\",\"clade11\",\"unknown\"),text.col=c(\"purple\",clade1col,clade2col,clade3col,clade4col,clade6col,clade8col,clade11col,\"black\"),bty=\"n\",xpd=T,inset = -0.05)
";

close RCMD;

print "Printed R commands to raxmlMLSThir/showTree$allname.R\n";


print SHOUT "
#cd raxmlMLSThir/
`R CMD BATCH --no-save --no-restore raxmlMLSThir/showTree$allname.R`
";


close SHOUT;

print "\nNext step:\nqsub run$allname

#./runMLSTanalysis2.pl -i odds1003dsts.xmfa
./runMLSTanalysis2.pl -i animaldsts.xmfa

module load phylip/3.69
getphyliptree.pl -i hirncycanimaldsts.mfa -n 1000

";

####################
### PAINT CHROMOSOMES SHOWING ALL COLORS
#####################


`mkdir paintNoNCYCprobe`;

open SHOUT, ">runfaChrompaint" or die "couldn't open runfaChrompaint : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runfaChrompaint
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=64gb

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3

cd  /lustre1/doudab/candidAnalysis/final/paintNoNCYCprobe
ln -s ../chr*.mfa .
rm chr1chr2*.mfa
rm chr*phy*
rm chr*_*

for strain in NCYC4144 NCYC4145 NCYC4146 NCYC597 12C 19F 1AA GC75 L26 P34048 P37005 P37037 P37039 P57055 P57072 P60002 P75010 P75016 P75063 P76055 P76067 P78042 P78048 P87 P94015 SC5314
do
  faChrompaint.pl -I chr -r \$strain -c ../clades -C ../colors -e 'NCYC 1AA SC5314_A'
done


";

close SHOUT;

print "\nNext step:\nqsub runfaChrompaint\n\n";


## FIGURE OUT 90% threshold for within clade divergence and visualise together with between clade divergence

open RCMD, ">paintNoNCYCprobe/divHist.R" or die "couldn't open paintNoNCYCprobe/divHist.R : $!";

print RCMD "

pdf('divHist.pdf')
options(scipen=999)

par(mfrow=c(2,1))
sameclade<-read.table(\"sameclade.pdiffs\")
head(sameclade)

quantile(sameclade\$V9,c(0.9,0.95,0.99))
summary(sameclade\$V9)

hist(sameclade\$V9,0:170/10000,col=\"blue\",xlab=\"Proportion of bases differing between pairs of sequences\",main=\"Within-clade comparisons\")
#abline(v=quantile(sameclade\$V9,0.95),col=\"red\")
abline(v=quantile(sameclade\$V9,0.9),col=\"green\")

diffclade<-read.table(\"diffclade.pdiffs\")
summary(diffclade\$V9)
head(diffclade)
quantile(diffclade\$V9,c(0.05,0.1,0.15,0.2))

hist(diffclade\$V9,0:170/10000,col=\"purple\",xlab=\"Proportion of bases differing between pairs of sequences\",main=\"Between-clade comparisons\")
abline(v=quantile(sameclade\$V9,0.9),col=\"green\")
#abline(v=quantile(sameclade\$V9,0.95),col=\"red\")
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

cd  /lustre1/doudab/candidAnalysis/final/paintNoNCYCprobe

cat *.pDiffs.tsv | awk '\$3==\$5{ print \$0; }' > sameclade.pdiffs  # only within clade diffs because \$3==\$5: oak strains and singletons excluded
cat *.pDiffs.tsv | grep -v pDiff | awk '\$5!=\$3 { print \$0; }' > diffclade.pdiffs

`R CMD BATCH --no-save --no-restore divHist.R`
";


close SHOUT;

print "\nNext step:\nqsub runDivhist\n\n";



########################################################
### PAINT CHROMOSOMES SHOWING ONLY COLORS with within clade similarity (<90% with -M option)
#######################################################


`mkdir paintNoNCYCprobeM90`;
#`mkdir paintWithNCYCprobeM90`;
`mkdir paintWithOakprobeM90`;

open SHOUT, ">runfaChrompaintM90" or die "couldn't open runfaChrompaintM90 : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runfaChrompaintM90
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=64gb

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3

#cd  /lustre1/doudab/candidAnalysis/final/paintNoNCYCprobeM90
#ln -s ../chr*.mfa .
#rm chr1chr2*.mfa
#rm chr*phy*
#rm chr*_*

#for strain in NCYC4144 NCYC4145 NCYC4146 NCYC597 12C 19F 1AA GC75 L26 P34048 P37005 P37037 P37039 P57055 P57072 P60002 P75010 P75016 P75063 P76055 P76067 P78042 P78048 P87 P94015 SC5314
#do
#  faChrompaint.pl -I chr -r \$strain -c ../clades -C ../colors -e 'NCYC 1AA SC5314_A'  -M $withinCladeMax
#done

#cd  /lustre1/doudab/candidAnalysis/final/paintWithNCYCprobeM90
cd  /lustre1/doudab/candidAnalysis/final/paintWithOakprobeM90
ln -s ../chr*.mfa .
rm chr1chr2*.mfa
rm chr*phy*
rm chr*_*

for strain in NCYC4144 NCYC4145 NCYC4146 NCYC597 12C 19F 1AA GC75 L26 P34048 P37005 P37037 P37039 P57055 P57072 P60002 P75010 P75016 P75063 P76055 P76067 P78042 P78048 P87 P94015 SC5314
do
  faChrompaint.pl -I chr -r \$strain -c ../clades -C ../colors -e 'NCYC597 1AA SC5314_A' -M $withinCladeMax
done


";

close SHOUT;

print "\nNext step:\nqsub runfaChrompaintM90\n\n";


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

cd  /lustre1/doudab/candidAnalysis/final/paintWithOakprobeM90

for strain in NCYC4144 NCYC4145 NCYC4146 NCYC597 12C 19F 1AA GC75 L26 P34048 P37005 P37037 P37039 P57055 P57072 P60002 P75010 P75016 P75063 P76055 P76067 P78042 P78048 P87 P94015 SC5314
do
  faChrompaint.pl -p \$strain.nearest.tsv -c ../clades -C ../colors2 
#  cat nearest\$strain.R | sed 's/mar=c(1,1,1,1)/mar=c(1,1,2.5,1)/g' | sed \"s/main=\\\"\$strain./main=\\\"/g\" | sed 's/.mfa\",yaxt=/\",cex.main=3,yaxt=/g' > newnearest\$strain.R
#  cat nearest\$strain.R | sed 's/mar=c(1,1,1,1)/mar=c(1,5,1,1),xpd=T/g' | sed \"s/main=\".*\",//g\" | sed \"s/xaxt='n')/xaxt='n',ylab='')\\ntext(-300000,50,labels=\\\"chr\\\",cex=3)/g\" > newnearest\$strain.R
  cat nearest\$strain.R | sed 's/mar=c(1,1,1,1)/mar=c(1,5.5,1,1),xpd=T/g' | sed \"s/main=\\\"\$strain./ylab='')\\ntext(-300000,50,labels=\\\"/g\" | sed \"s/.mfa\\\")/\\\",cex=3.5)/g\" > newnearest\$strain.R
 


  R CMD BATCH --no-restore newnearest\$strain.R
  mv \$strain.nearest.pdf \${strain}"."newnearest.pdf
done


";

close SHOUT;

print "\nNext step:\nqsub runRePaint\n\n";



#########################
### DRAW PHYLOGENETIC FIGURE
##########################

# generate phylogenies for small genomic regions (code adapted from runPhyloSubseq.pl

my @chrs = qw(R 5);
my @starts = qw(2000000 150000);
my @ends = qw(2100000 250000);

open SHOUT, ">runSubPhylo" or die "couldn't open runSubPhylo : $!";

printphyloheader("runSubPhylo");

for (my $i=0; $i < @chrs; $i++) {
	
	my $outname = "chr$chrs[$i]"."_$starts[$i]"."_$ends[$i]";

	`mkdir raxml$outname`;
	`grep ">" chr$chrs[$i].mfa | sed 's/>//g' > names`;

## CREATE ALIGNMENT IN PHYLIP FORMAT FOR EACH CHR WITH LOWERCASE LETTERS CONVERTED TO Ns and Ns added to the end of chr (with alcat.pl)
## THEN RUN RAXML

	print SHOUT"
faChooseSubseq.pl -i chr$chrs[$i].mfa -n names -s '$starts[$i]..$ends[$i]' -o $outname.mfa 
cat $outname.mfa | awk -F: '{ print \$1; }' > temp
mv temp $outname.mfa
fa2phylip.pl -i $outname.mfa > $outname.phy


rm raxml$outname/*chr*
raxmlHPC-PTHREADS-AVX -T 4 -f a -x 12345 -p 12345 -m GTRGAMMA -N $bootstraps -s $outname.phy -n $outname -w $dir/raxml$outname/

";				

	print "\ncreated commands to generate an alignment and raxml phylogeny for $chrs[$i] $starts[$i]..$ends[$i] in $outname.mfa\n"; 


}
close SHOUT;

print "\nNext step:\nqsub runSubPhylo\n\n";

## DRAW 4-PART PHYLOGENETIC FIGURE : genome-wide, chr3, subphylogenies


open RCMD, ">phyloFig.R" or die "couldn't open phyloFig.R : $!";

print RCMD "

library(ape)
pdf('phyloFig.pdf')

par(mfrow=c(2,2),mar=c(1,1,1,1))
";

printRcmds4fig("chr1chr2chr3chr4chr5chr6chr7chrR");

print RCMD "
legend(\"bottomleft\", \"a.\", bty=\"n\")
legend(\"bottomright\",c(\"oak\",\"clade1\",\"clade2\",\"clade3\",\"clade4\",\"clade6\",\"clade8\",\"clade11\"),text.col=c(\"purple\",\"#dd1c77\",\"orange\",\"blue\",\"#bdbdbd\",\"#fa9fb5\",\"dark grey\",\"#5aae61\"),bty=\"n\",cex=0.8)
";

printRcmds4fig("chr3");

print RCMD "legend(\"bottomleft\", \"b.\", bty=\"n\")\n";

my @legend = qw(c d);

for (my $i=0; $i < @chrs; $i++) {
	
	my $outname = "chr$chrs[$i]"."_$starts[$i]"."_$ends[$i]";

	printRcmds4fig($outname);

	print RCMD "legend(\"bottomleft\", \"$legend[$i].\", bty=\"n\")\n";
}
close RCMD;

print "\nPrinted R commands for drawing 4-part Figure that shows genome-wide phylogeny and incongruence to phyloFig.R\n";

open SHOUT, ">runPhyloFig" or die "couldn't open runPhyloFig : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runPhyloFig
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=64gb

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3

cd /lustre1/doudab/candidAnalysis/final/

`R CMD BATCH --no-save --no-restore phyloFig.R`

";

close SHOUT;

print "\nNext step to run the R commands:\nqsub runPhyloFig\n\n";



## DRAW 2-PART PHYLOGENETIC FIGURE : genome mlst, animal mlst with RAXML

open RCMD, ">phyloFig2part.R" or die "couldn't open phyloFig2part.R : $!";

print RCMD "library(ape)\n";

if ($tiff eq "no") { print RCMD "pdf('phyloFig2part.pdf',height = 7, width = 5)\n"; }
else { print RCMD "bitmap(\"phyloFig2part.tiff\", height = 7, width = 4, units = 'in', type = \"tiff24nc\", res = 300)\n"; }

print RCMD "

par(mfrow=c(2,1),cex=1.2)
";

printRcmds4fig("MLSThir");

print RCMD "
par(mar=c(3,3,1,3),cex=1,xpd=T)
legend(\"bottomleft\", \"a.\", bty=\"n\",cex=2,inset = -0.2)
legend(\"bottomright\",c(\"oak\",\"clade1\",\"clade2\",\"clade3\",\"clade4\",\"clade6\",\"clade8\",\"clade9\",\"clade10\",\"clade11\",\"unknown\"),text.col=c(wildcol,clade1col,clade2col,clade3col,clade4col,clade6col,clade8col,clade9col,clade10col,clade11col,ncyccol),bty=\"n\",cex=1.5,xpd=T,inset = -0.1)
";

printRcmds4fig("hirncycanimaldsts");

print RCMD "
par(mar=c(3,3,1,1),cex=1,xpd=T)
legend(\"bottomleft\", \"b.\", bty=\"n\",cex=2, inset = -0.1, xpd=T)
";

close RCMD;

print "\nPrinted R commands for drawing 2-part Figure that shows genome-wide phylogeny and one other phyloFig2part.R\n";

open SHOUT, ">runPhyloFig2part" or die "couldn't open runPhyloFig2part : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runPhyloFig2part
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=64gb

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3

cd /lustre1/doudab/candidAnalysis/final/

`R CMD BATCH --no-save --no-restore phyloFig2part.R`

";

if ($tiff eq "yes") { print SHOUT "`tiff2pdf phyloFig2part.tiff -o phyloFig2part.pdf`\n"; }

close SHOUT;

print "\nNext step to run the R commands:\nqsub runPhyloFig2part\n\n";

## DRAW 1-PART animal FIGURE :  animal mlst with RAXML

open RCMD, ">phyloFigAnimal.R" or die "couldn't open phyloFigAnimal.R : $!";

print RCMD "library(ape)\n";

if ($tiff eq "no") { print RCMD "pdf('phyloFigAnimal.pdf',height = 7, width = 5)\n"; }
else { print RCMD "bitmap(\"phyloFigAnimal.tiff\", height = 7, width = 4, units = 'in', type = \"tiff24nc\", res = 300)\n"; }

print RCMD "

par(mfrow=c(1,1))
";

printRcmds4fig("hirncycanimaldsts");

print RCMD "
par(mar=c(3,3,1,1),cex=1,xpd=T)
legend(\"bottomleft\",c(\"oak\",\"clade1\",\"clade2\",\"clade3\",\"clade4\",\"clade6\",\"clade8\",\"clade9\",\"clade10\",\"clade11\",\"unknown\"),text.col=c(wildcol,clade1col,clade2col,clade3col,clade4col,clade6col,clade8col,clade9col,clade10col,clade11col,ncyccol),bty=\"n\",cex=1.5,xpd=T)
";


close RCMD;

print "\nPrinted R commands for drawing 1-part Figure that shows animal phylogeny and one other phyloFigAnimal.R\n";

open SHOUT, ">runPhyloFigAnimal" or die "couldn't open runPhyloFigAnimal : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runPhyloFigAnimal
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=64gb

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3

cd /lustre1/doudab/candidAnalysis/final/

`R CMD BATCH --no-save --no-restore phyloFigAnimal.R`

";

if ($tiff eq "yes") { print SHOUT "`tiff2pdf phyloFigAnimal.tiff -o phyloFigAnimal.pdf`\n"; }

close SHOUT;

print "\nNext step to run the R commands:\nqsub runPhyloFigAnimal\n\n";


## DRAW 1-PART PHYLOGENETIC FIGURE : genome-wide

open RCMD, ">phyloFig1part.R" or die "couldn't open phyloFig1part.R : $!";

print RCMD "library(ape)\n";

if ($tiff eq "no") { print RCMD "pdf('phyloFig1part.pdf')\n"; }
else { print RCMD "bitmap(\"phyloFig1part.tiff\", type = \"tiff24nc\", res = 300)\n"; }

print RCMD "

par(mfrow=c(1,1),mar=c(3,3,3,3),cex=1,xpd=T)
";

printRcmds4fig("chr1chr2chr3chr4chr5chr6chr7chrR");

print RCMD "
par(mar=c(1,1,1,1),cex=1)
#legend(\"bottomleft\", \"a.\", bty=\"n\",cex=2)
legend(\"bottomright\",c(\"oak\",\"clade1\",\"clade2\",\"clade3\",\"clade4\",\"clade6\",\"clade8\",\"clade11\",\"unknown\"),text.col=c(wildcol,clade1col,clade2col,clade3col,clade4col,clade6col,clade8col,clade11col,\"brown\"),bty=\"n\",cex=1.2,inset = -0.05)
";

close RCMD;

print "\nPrinted R commands for drawing 1-part Figure that shows genome-wide phylogeny only\n";

open SHOUT, ">runPhyloFig1part" or die "couldn't open runPhyloFig1part : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runPhyloFig1part
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=64gb

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3

cd /lustre1/doudab/candidAnalysis/final/

`R CMD BATCH --no-save --no-restore phyloFig1part.R`

";

if ($tiff eq "yes") { print SHOUT "`tiff2pdf phyloFig1part.tiff -o phyloFig1part.pdf`\n"; }

close SHOUT;

print "\nNext step to run the R commands:\nqsub runPhyloFig1part\n\n";



## DRAW 2-PART INCONGRUENCE FIGURE : chr3, subphylogeny


open RCMD, ">phyloFig2.R" or die "couldn't open phyloFig2.R : $!";

print RCMD "

library(ape)
";

if ($tiff eq "no") { print RCMD "pdf('phyloFig2.pdf',height = 4, width = 7,paper=\"special\")\n"; }
#else { print RCMD "bitmap(\"phyloFig2.tiff\", height = 4, width = 7, units = 'in', type = \"tiff24nc\", res = 300)\n"; }
else { print RCMD "bitmap(\"phyloFig2.tiff\",type = \"tiff24nc\", res = 300)\n"; }


print RCMD "
#par(mfrow=c(1,2),mar=c(1,1,1,1),cex=1.2)
par(mfrow=c(1,1),mar=c(1,1,1,1),cex=1)

";

#printRcmds4fig("chr7");

print RCMD "
#legend(\"bottomleft\", \"a.\", bty=\"n\",cex=1.5)
";


my $outname = "chr5_150000_250000";
printRcmds4fig($outname);

print RCMD "
#legend(\"bottomleft\", \"b.\", bty=\"n\",cex=1.5)
legend(\"topright\",c(\"oak\",\"clade1\",\"clade2\",\"clade3\",\"clade4\",\"clade6\",\"clade8\",\"clade11\",\"unknown\"),text.col=c(wildcol,clade1col,clade2col,clade3col,clade4col,clade6col,clade8col,clade11col,ncyccol),bty=\"n\",cex=1.5)
";

close RCMD;

#print "\nPrinted R commands for drawing 2-part Figure that shows incongruence to phyloFig2.R\n";
print "\nPrinted R commands for drawing 1-part Figure that shows incongruence to phyloFig2.R (historic name)\n";


open SHOUT, ">runPhyloFig2" or die "couldn't open runPhyloFig2 : $!";

print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N runPhyloFig2
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=64gb

#PBS -M doudab\@uga.edu 
#PBS -m ae

module load R/3.2.3

cd /lustre1/doudab/candidAnalysis/final/

`R CMD BATCH --no-save --no-restore phyloFig2.R`

";

if ($tiff eq "yes") { print SHOUT "`tiff2pdf phyloFig2.tiff -o phyloFig2.pdf`\n"; }

close SHOUT;

print "\nNext step to run the R commands:\nqsub runPhyloFig2\n\n";





#########################
### DRAW HETEROZYGOSITY FIGURE 2
##########################

# Print Rcmds to run with R CMD BATCH
open RCMD, ">$hethist.R" or die "couldn't open $hethist.R : $!";

if ($tiff eq "no") { print RCMD "pdf(\"$hethist.pdf\", height = 3, width = 7)\n"; }
else { print RCMD "bitmap(\"hethist.tiff\", height = 2, width = 7, units = 'in', type = \"tiff24nc\", res = 300)\n"; }

print RCMD "
par(mfrow=c(1,3),mar=c(5, 4, 2, 2),cex=1.5);
";

#foreach my $strain qw(NCYC4144 NCYC4145 NCYC4146 NCYC597 SC5314 1AA) {
foreach my $strain qw(NCYC4144 NCYC4145 NCYC4146) {


	print RCMD "

rm(list=ls())
#data<-read.table(\"$strain.I.calls\",header=T)
#attach(data)
#head(data)
#data2<-read.table(\"$strain.I.vcf.all\",header=T)
#head(data2)
";

	# READ IN DATA FOR THE HISTOGRAM FROM THE FIRST RUN OF vcf2alleles.pl

	my ($gettingclose, @glhet_snp, @glhet_length);
	open IN, "<$strain.I.Rcmds.Rout" or die "$strain.I.Rcmds.Rout : $!";
	while (<IN>) {
		if (/cbind\(glhet_snp,glhet_length\)/) { $gettingclose++; }
		if ((defined $gettingclose) && (/^\s*\[[0-9]+,\]\s+(\d+)\s+(\d+)\s*$/m)) { 
		#	print "$strain.I.Rcmds.Rout\tglhet_snp: $1, glhet_length: $2\n"; 
			push (@glhet_snp,$1);
			push (@glhet_length,$2);		
		}
	}
	close IN;

	print RCMD "glhet_snp <- c(".join( ',', @glhet_snp ).")\n";
	print RCMD "glhet_length <- c(".join( ',', @glhet_length ).")\n";

	print RCMD "hist(glhet_snp/glhet_length,breaks=0:32/1000,xlab=\"$strain heterozygosity\",col=\"blue\",main=\"\")\n";
	print RCMD "abline(v=$maxLOH,col=\"red\")\n";
}

#########################
### DRAW HETEROZYGOSITY SUPP FIGURE 
##########################


print RCMD "
pdf(\"$hethist"."Supp.pdf\")
par(mfrow=c(4,3));
";

print "prepare to draw $hethist for the following strains:\n"; 
foreach my $strain (sort @strains) {

	if (($strain =~ /NCYC4/)) { next; }	# skip strains that were in the main figure
#	if (($strain =~ /NCYC/) || ($strain eq "SC5314") || ($strain eq "1AA")) { next; }	# skip strains that were in the main figure
	print "$strain ";
	print RCMD "

rm(list=ls())
#data<-read.table(\"$strain.I.calls\",header=T)
#attach(data)
#head(data)
#data2<-read.table(\"$strain.I.vcf.all\",header=T)
#head(data2)
";

	# READ IN DATA FOR THE HISTOGRAM FROM THE FIRST RUN OF vcf2alleles.pl

	my ($gettingclose, @glhet_snp, @glhet_length);
	open IN, "<$strain.I.Rcmds.Rout" or die "$strain.I.Rcmds.Rout : $!";
	while (<IN>) {
		if (/cbind\(glhet_snp,glhet_length\)/) { $gettingclose++; }
		if ((defined $gettingclose) && (/^\s*\[[0-9]+,\]\s+(\d+)\s+(\d+)\s*$/m)) { 
		#	print "$strain.I.Rcmds.Rout\tglhet_snp: $1, glhet_length: $2\n"; 
			push (@glhet_snp,$1);
			push (@glhet_length,$2);		
		}
	}
	close IN;

	print RCMD "glhet_snp <- c(".join( ',', @glhet_snp ).")\n";
	print RCMD "glhet_length <- c(".join( ',', @glhet_length ).")\n";

	print RCMD "hist(glhet_snp/glhet_length,breaks=0:35/1000,xlab=\"$strain heterozygosity\",col=\"blue\",main=\"\")\n";
	print RCMD "abline(v=$maxLOH,col=\"red\")\n";
}
print "\n";

close RCMD;

# RUN THE R COMMANDS WITH OUTPUTS GOING OUT TO THE DEFAULT FILE NAMES

open SHOUT, ">runHethist" or die "couldn't open runHethist : $!";

printheader("runHethist");

print SHOUT "`R CMD BATCH --no-save --no-restore $hethist.R`\n";

if ($tiff eq "yes") { print SHOUT "`tiff2pdf hethist.tiff -o hethist.pdf`\n"; }

close SHOUT;

print "\nNext step to run the R commands:\nqsub runHethist\n";

print "R output will then be in in $hethist.Rout and histograms of heterozygosity will be in $hethist.pdf\n\n";


###########################################################
## CODE ADAPTED FROM runHetanalysis.pl for TABLE 2 AND HETEROZYGOSITY ANALYSES
###########################################################

open TABLE, ">$table" or die "couldn't open $table: $!";

opendir (DIR, ".") or die "couldn't open . : $!";

#my (%strains,@hetfiles);					# collect the names of het summary files and strains 
my (@hetfiles);					# collect the names of het summary files and strains 

while (defined (my $name = readdir(DIR))) {

	if ($name =~ /^(\S+)\.simple2d\.het\.txt$/i) { 
		print "$name\tstrain:$1\n"; 	
		$strains{$name} = $1;
		push (@hetfiles, $name);
	}

}		

# EXTRACT GENOMEWIDE HETEROZYGOSITY DATA FROM *\.simple2d\.het\.txt files


open DATA, ">$datafile" or die "couldn't open $datafile : $!";
print DATA "strain\tfile\tmtl\tclade\tq40ps\tq40l\tq40het\tfps\tfl\tfhet\tfdps\tfdl\tfdhet\tlohl\n";
print TABLE "Strain & MTL & Clade & Heterozygosity\\footnote\\ & LOH length (Mbp)\\footnote\\ & Filtered length (Mbp)\\footnote\\ & Filtered heterozygous sites & Filtered heterozygosity\\footnote\\ \\\\ \\hline\n";


my (%hqh2strains,%hqh,%lohl,%dfl,%dfh,%dfps,%hqps,%hql,%fl,%fps,%fh);

my $nncyc = 0;
my $sumhqhncyc = 0;
my $sumhqhclinic = 0;
my $nclinic = 0;
my $sumlohncyc = 0;
my $sumlohclinic = 0;
my $sumdfhncyc = 0;
my $sumdfhclinic = 0;
foreach my $file (@hetfiles) { 

	unless (defined $mtl{$strains{$file}}) { $mtl{$strains{$file}} = "?"; }
	unless (exists $clades{$strains{$file}}) { $clades{$strains{$file}} = "?"; }

	print DATA "$strains{$file}\t$file\t$mtl{$strains{$file}}\t$clades{$strains{$file}}";

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
	print "extracted LOH data for $strains{$file} from $gffile.\n";
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

my $nstrains = scalar keys %hqh2strains;
foreach my $het (reverse sort keys %hqh2strains) {
	my $strain = $hqh2strains{$het};

	my $mtl;
	if ($mtl{$strain} eq "het") { $mtl = "\$a/\\alpha\$"; }
	elsif ($mtl{$strain} eq "homa") { $mtl = "\$a/a\$"; }
	elsif ($mtl{$strain} eq "homalpha") { $mtl = "\$\\alpha/\\alpha\$"; }
	else { $mtl = $mtl{$strain}; }

	print TABLE "$strain & $mtl & $clades{$strain} & ".sprintf("%.4f",$hqh{$strain})." & ".sprintf("%.1f",$lohl{$strain}/1000000)." & ".sprintf("%.1f",$dfl{$strain}/1000000)." & $dfps{$strain} & ".sprintf("%.4f",$dfh{$strain})." \\\\ ";

	if ($strain eq "NCYC4144") { 
		print TABLE "\n\\textbf{3 oak strains} &  & Mean: & ".sprintf("%.4f",$sumhqhncyc/$nncyc)."*** & ".sprintf("%.1f",($sumlohncyc/$nncyc)/1000000)."* & & &".sprintf("%.4f",$sumdfhncyc/$nncyc)."** \\\\ \\\\ \n";
	}
	elsif ($strain eq "1AA") {			# is this the last entry?
		print TABLE "\n\\textbf{".($nstrains-4)." clinical strains}\\footnote\\ &  & Mean: & ".sprintf("%.4f",$sumhqhclinic/$nclinic)."*** & ".sprintf("%.1f",($sumlohclinic/$nclinic)/1000000)."* & & & ".sprintf("%.4f",$sumdfhclinic/$nclinic)."** \\\\ \\hline \n";
	}
	else { print TABLE "\n"; }


}
close TABLE;

print "\nPrinted a data for all $nstrains strains from ".@hetfiles." *.simple2d.het.txt files to $table.\n\n";


open RCMD, ">Hetanalysis.Rcmd" or die "couldn't open Hetanalysis.Rcmd : $!";

print RCMD "

rm(list=ls())
data<-read.table(\"$datafile\",header=T)
data2<-data[data\$strain!=\"1AA\",]	# drop 1AA because it is not independent of SC5314
attach(data2)
head(data2)

oakq40het<-q40het[grep(\"NCYC4\",strain)]
clinicq40het<-q40het[grep(\"NCYC4\",strain,invert=T)]		# this includes NCYC597 and WT_SC5314

oakq40l<-q40l[grep(\"NCYC4\",strain)]
clinicq40l<-q40l[grep(\"NCYC4\",strain,invert=T)]
oakq40l
sort(clinicq40l)
summary(oakq40l)
summary(clinicq40l)
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

summary(fdl)
sort(fdl)
length(oakfdhet)
length(clinicfdhet)
summary(oakfdhet)
summary(clinicfdhet)
data2[grep(\"NCYC4\",strain),]
wilcox.test(oakfdhet,clinicfdhet)\

o<-order(fdhet)
strain[o]
fdhet[o]

# OAK HET STRAINS VS CLINICAL HET STRAINS

hetdata2<-data2[mtl==\"het\",]
oakfdhethet<-hetdata2\$fdhet[grep(\"NCYC4\",hetdata2\$strain)]
clinicfdhethet<-hetdata2\$fdhet[grep(\"NCYC4\",hetdata2\$strain,invert=T)]		# this includes NCYC597 and WT_SC5314

length(oakfdhethet)
length(clinicfdhethet)
summary(oakfdhethet)
summary(clinicfdhethet)
wilcox.test(oakfdhethet,clinicfdhethet)

# CLINICAL HET STRAINS VS CLINICAL HOM STRAINS

homdata2<-data2[mtl!=\"het\",]
head(homdata2)
clinicfdhethom<-homdata2\$fdhet[grep(\"NCYC4\",homdata2\$strain,invert=T)]		# this includes NCYC597 and WT_SC5314

length(clinicfdhethet)
length(clinicfdhethom)
summary(clinicfdhethet)
summary(clinicfdhethom)
wilcox.test(clinicfdhethet,clinicfdhethom)


# DO OAK STRAINS HAVE LESS LOH?						

oaklohl<-lohl[grep(\"NCYC4\",strain)]
cliniclohl<-lohl[grep(\"NCYC4\",strain,invert=T)]		# this includes NCYC597 and WT_SC5314

length(oaklohl)
length(cliniclohl)
summary(oaklohl)
summary(cliniclohl)
data2[grep(\"NCYC4\",strain),]
wilcox.test(oaklohl,cliniclohl)
";


close RCMD;

open SHOUT, ">runHetanalysis" or die "couldn't open runHetanalysis : $!";

printheader("runHetanalysis");

print SHOUT "`R CMD BATCH --no-save --no-restore Hetanalysis.Rcmd`\n";

close SHOUT;

print "\nNext step to run the R commands:\nqsub runHetanalysis\n";
print "Then results of analysis will be in in Hetanalysis.Rcmd.Rout, $datafile and $table\n\n";


print "NB to see similarity with ST sequences matching oak:\ncat hirncycanimaldsts.mfa SToakseqs.fas > STsmatchingOak.mfa\n\n";


##
## DAVE HALL's SUGGESTION : heterozygosity analysis in zero LOH regions
##

## REVIEW and pool LOH regions for all strains

my (%lohstarts,%lohends);

foreach my $file (@callfiles) {						# read annotations for all files

	my $gffile = "$strains{$file}.LOH.I.gff";

	if ($strains{$file} eq "1AA") { next; }

	if (defined $gffile) { 			# Read annotations if they exist
		foreach my $chr (@chr) {
			if (defined $atype{"$file$chr"}[0]) {
				my @atype = @{$atype{"$file$chr"}};
				my @astart = @{$astart{"$file$chr"}};
				my @aend = @{$aend{"$file$chr"}};

	#			print "$file\t$gffile\t$chr :\n";
				for (my $i=0; $i<@astart; $i++) {	# save all the LOH descriptions for every chromosome
					push(@{$lohstarts{$chr}}, $astart[$i]); 
					push(@{$lohends{$chr}}, $aend[$i]); 
				}

			}		
		}
	}
}

my (%mergedstarts,%mergedends);

foreach my $chr (@chr) {
	my @lohstarts = @{$lohstarts{$chr}};
	my @lohends = @{$lohends{$chr}};
	my %seenstart;

	my $inasegment = "no";

	print "$chr\t";


	for (my $i=0; $i < @lohstarts; $i++) {
		if (!defined $mergedstarts{$chr}) {	# the first LOH region for this chromosome
			push(@{$mergedstarts{$chr}}, $lohstarts{$chr}[$i]); 
			push(@{$mergedends{$chr}}, $lohends{$chr}[$i]); 
			next;
		}
		for (my $j=0; $j < @{$mergedstarts{$chr}}; $j++) {

							# contained in a segment
			if ( ($lohstarts{$chr}[$i] > $mergedstarts{$chr}[$j]) && ($lohends{$chr}[$i] < $mergedends{$chr}[$j]) ) { 
				$inasegment = "yes";
			}
							# need to extend the LOH start for this segment
			elsif ( ($lohstarts{$chr}[$i] < $mergedstarts{$chr}[$j]) && ($lohends{$chr}[$i] > $mergedstarts{$chr}[$j]) ) { 
				$mergedstarts{$chr}[$j] = $lohstarts{$chr}[$i];
			#	print "\tnew start: $mergedstarts{$chr}[$j]..$mergedends{$chr}[$j]\n";
				$inasegment = "yes";
 			}

							# need to extend the LOH end for this segment
			elsif ( ($lohends{$chr}[$i] > $mergedends{$chr}[$j]) && ($lohstarts{$chr}[$i] < $mergedends{$chr}[$j]) ) { 
				$mergedends{$chr}[$j] = $lohends{$chr}[$i];
			#	print "\tnew end: $mergedstarts{$chr}[$j]..$mergedends{$chr}[$j]\n";
				$inasegment = "yes";
			}

							# need to extend the LOH start AND end for this segment
			elsif ( ($lohends{$chr}[$i] > $mergedends{$chr}[$j]) && ($lohstarts{$chr}[$i] < $mergedstarts{$chr}[$j]) ) { 
				$mergedstarts{$chr}[$j] = $lohstarts{$chr}[$i];
				$mergedends{$chr}[$j] = $lohends{$chr}[$i];
				$inasegment = "yes";
			}

						
		}
		if ($inasegment eq "yes") { $inasegment = "no"; next; }
		elsif ($inasegment eq "no") {
			unless ($seenstart{$lohstarts{$chr}[$i]}++) { 
				push(@{$mergedstarts{$chr}}, $lohstarts{$chr}[$i]);
				push(@{$mergedends{$chr}}, $lohends{$chr}[$i]);
			#	print "\tnew segment:  $lohstarts{$chr}[$i]..$lohends{$chr}[$i] (starts: @{$mergedstarts{$chr}}\n";
			}
		}		
					
	}

	print "$chr after merger: ";
	for (my $j=0; $j < @{$mergedstarts{$chr}}; $j++) {
		print "$mergedstarts{$chr}[$j]..$mergedends{$chr}[$j] ";
	}
	print "\n";
}

# manually summarise non-LOH regions (approx 1.2 Mbp total on in 5 regions on 4 chromosomes)

my %nonLOHchr = (
	1 => "chr1",	
	2 => "chr2",	
	3 => "chr2",
	4 => "chr4",	
	5 => "chr6",		
);

my %nonLOHstart = (
	1 => 1500000,	
	2 => 900000,
	3 => 1700000,
	4 => 800000,
	5 => 900000,					
);

my %nonLOHend = (
	1 => 1700000,			
	2 => 1000000,	
	3 => 2100000,	
	4 => 1200000,
	5 => 1000000,		
);

#####################################################################
## COMPARE HETEROZYGOSITY IN REGIONS THAT ARE NON-LOH FOR ALL STRAINS
#####################################################################

open RCMD, ">HetanalysisDH.Rcmd" or die "couldn't open HetanalysisDH.Rcmd : $!";

my $k=0;

print RCMD "

rm(list=ls())
strains<-\"none\"
totalhet<-0
totallen<-0
";
foreach my $strain (sort @strains) {
#foreach my $strain ("12C") {

	print "$strain ";


	print RCMD "
data<-read.table(\"$strain.simple2d.calls\",header=T)
attach(data)
head(data)

sum(pALT<=0.8&pALT>=0.2&QUAL>40)

";

	print RCMD "

data2<-read.table(\"$strain.simple2d.all\",header=T)
head(data2)

table(data\$filter[QUAL>40])
table(data2\$filter[data2\$QUAL>=40])
table(data2\$filter)

table(data2\$filter[data2\$QUAL>40]) 	# these are the data we want!

het<-0
len<-0
";


	for (my $i=1; $i<=5; $i++) {
		print RCMD "
het[$i]<-sum(pALT<=0.8&pALT>=0.2&QUAL>40&chr==\"$nonLOHchr{$i}\"&pos>$nonLOHstart{$i}&pos<=$nonLOHend{$i})
len[$i]<-sum(data2\$QUAL>40&data2\$chr==\"$nonLOHchr{$i}\"&data2\$pos>$nonLOHstart{$i}&data2\$pos<=$nonLOHend{$i}&data2\$filter==\"no\")
";
	}

	print RCMD "
# summary for $strain
het
len
het/len
sum(het)/sum(len)
";
	if ($strain ne "1AA") {		# exclude 1AA from the wilcox.test
		$k++;
		print RCMD "	
strains[$k]<-\"$strain\"
totalhet[$k]<-sum(het)
totallen[$k]<-sum(len)

strains
totalhet
totallen
";
	}


}


print RCMD "

# SUMMARISE POOLED DATA FOR THESE 5 LOCI WITH WILCOX TEST BETWEEN OAK AND NON_OAK (excluding 1AA)

strains
totalhet
totallen

oakhet<-totalhet[grep(\"NCYC4\",strains)]
oaklen<-totallen[grep(\"NCYC4\",strains)]


nonoakhet<-totalhet[grep(\"NCYC4\",strains,invert=T)]
nonoaklen<-totallen[grep(\"NCYC4\",strains,invert=T)]

sort(oakhet/oaklen)
sort(nonoakhet/nonoaklen)

summary(oakhet/oaklen)
summary(nonoakhet/nonoaklen)

wilcox.test(oakhet/oaklen,nonoakhet/nonoaklen)

sort(oaklen)
sort(nonoaklen)
wilcox.test(oaklen,nonoaklen)

";


close RCMD;

open SHOUT, ">runHetanalysisDH" or die "couldn't open runHetanalysisDH : $!";

printheader("runHetanalysisDH");

print SHOUT "`R CMD BATCH --no-save --no-restore HetanalysisDH.Rcmd`\n";

close SHOUT;

print "\nNext step to run the R commands:\nqsub runHetanalysisDH\n";
print "Then results of analysis will be in in HetanalysisDH.Rcmd.Rout\n\n";



#####################################################################
## COMPARE HETEROZYGOSITY IN REGIONS THAT ARE *** UNFILTERED *** FOR ALL STRAINS
#####################################################################


# SUMMARISE ALL THE POSITIONS THAT NEED TO BE FILTERED

my (%allfpos,%seenallf,%chrmax);

print "Reading *.simple2d.all into a combined filter ..\n";
foreach my $strain (sort @strains) {
#foreach my $strain qw(12C 19F) {

	next;	# for a fast run (for a full run, comment this out)

	print "$strain \n";

	if ($strain eq "1AA") {	print "omitting $strain from filter\n"; next; }	# exclude 1AA from the filter



	open ALLSITES, "<$strain.simple2d.all" or die "couldn't open $strain.simple2d.all : $%!";
	while (<ALLSITES>){

		if (/(chr\S+)\s+(\d+)\s+(\S+)\s+(\S+)/) { 
			my $chr = $1;
			my $pos = $2,
			my $qual = $3;
			my $filter = $4;

			if ( (!defined $chrmax{$chr}) || ($pos > $chrmax{$chr}) ) { $chrmax{$chr} = $pos; } 
			if (($filter eq 'no') && ($qual > 40)) { next; }	
			
			unless ($seenallf{$chr."_$pos"}++) { push(@{$allfpos{$chr}},$pos); }		
		}
	}

}

my $allfilteredpos = keys %seenallf;

my $totallength = 0;

foreach my $chr (@chr) {
	print "$chr\t$chrmax{$chr} bp\n";
	$totallength += $chrmax{$chr};
}
print "\nTotal length: $totallength bp\n";
print "Unfiltered length: ".($totallength-$allfilteredpos)." bp\n";

print "\nNumber of filtered positions: ".$allfilteredpos."\n\nReading and printing a summary of heterozygosity for every strain one chromosome at a time ..\n";

foreach my $strain (sort @strains) {
#foreach my $strain qw(12C 19F) {

	next;	# for a fast run (for a full run, comment this out)


	open OUT, ">$strain.allfhet.txt" or die "couldn't open $strain.allfhet.txt : $!";
	print OUT "chr\tpos\tpALT\tfiltered\n";

	foreach my $chr (@chr) {
		print "$chr\t$chrmax{$chr} bp\n";
		$totallength += $chrmax{$chr};

		my (%fpos4thischr,%pAlt4thischr);
		$fpos4thischr{$_}++ for (@{$allfpos{$chr}});

		open HETSITES, "<$strain.simple2d.calls";

		while (<HETSITES>) {
			if (/$chr\s+(\d+)\s+\S+\s+\S+\s+\S+\s+\d+\s+\d+\s+\d+\s+\d+\s+(\S+)/) { $pAlt4thischr{$1} = $2; }
		}

		for (my $i=1; $i <= $chrmax{$chr}; $i++) {

			if (!defined $pAlt4thischr{$i}) { $pAlt4thischr{$i} = 0; }
			if (!defined $fpos4thischr{$i}) { $fpos4thischr{$i} = "no"; }
			else { $fpos4thischr{$i} = "yes"; }

			print OUT "$chr\t$i\t$pAlt4thischr{$i}\t$fpos4thischr{$i}\n";
		}
	}
	close OUT;
	print "printed summary is in $strain.allfhet.txt\n";
}


open RCMD, ">HetAllUnfilt.Rcmd" or die "couldn't open HetAllUnfilt.Rcmd : $!";

print RCMD "

rm(list=ls())
strains<-\"none\"
het<-0
len<-0
";

my $i = 0;

foreach my $strain (sort @strains) {
#foreach my $strain qw(12C 19F) {


	print "$strain ";


	print RCMD "
data<-read.table(\"$strain.allfhet.txt\",header=T)
attach(data)
head(data)

table(filtered)
sum(pALT<=0.8&pALT>=0.2&filtered==\"no\")

";

	if ($strain eq "1AA") { next; }	

	$i++;

	print RCMD "
strains[$i]<-\"$strain\"
het[$i]<-sum(pALT<=0.8&pALT>=0.2&filtered==\"no\")
len[$i]<-sum(filtered==\"no\")

strains
het
het/len
";

}


print RCMD "

# SUMMARISE POOLED DATA FOR ALL STRAINS WITH WILCOX TEST BETWEEN OAK AND NON_OAK (excluding 1AA)

strains
het
len

oakhet<-het[grep(\"NCYC4\",strains)]
oaklen<-len[grep(\"NCYC4\",strains)]


nonoakhet<-het[grep(\"NCYC4\",strains,invert=T)]
nonoaklen<-len[grep(\"NCYC4\",strains,invert=T)]

sort(oakhet/oaklen)
sort(nonoakhet/nonoaklen)

summary(oakhet/oaklen)
summary(nonoakhet/nonoaklen)

wilcox.test(oakhet/oaklen,nonoakhet/nonoaklen)

";


close RCMD;

open SHOUT, ">runHetAllUnfilt" or die "couldn't open runHetAllUnfilt : $!";

printheader("runHetAllUnfilt");

print SHOUT "`R CMD BATCH --no-save --no-restore HetAllUnfilt.Rcmd`\n";

close SHOUT;

print "\nNext step to run the R commands:\nqsub runHetAllUnfilt\n";
print "Then results of analysis will be in in HetAllUnfilt.Rcmd.Rout\n\n";


#####################
## REMAKE TABLE 2 to include Dave Hall's suggested analysis
## and make Table S1 to make room for an expansion of Table 2
######################

open IN, "<HetAllUnfilt.Rcmd.Rout" or die "couldn't open HetAllUnfilt.Rcmd.Rout : $!";

my $Rin;
while (<IN>) { $Rin .= $_; }	# from the Routput

my (%allhps,$allength);		# read the length of unfiltered sequence ($allength) and the number of heterozygous sites for each strain (%allhps)
while ($Rin =~ /data<-read.table\(\"(\S+?)\.allfhet\.txt.+?table\(filtered\).*?(\d+)\s+\d+.*?\[1\]\s+(\d+)/gsc) {
	$allength = $2;
	$allhps{$1} = $3;
	print "$1\tlength: $2\thetps: $3\thet: ".($3/$2)."\n";
}

open TABLE, ">$tableDH" or die "couldn't open $tableDH: $!";
open STABLE, ">$stable" or die "couldn't open $stable : $!";

print STABLE "Strain\tMTL\tClade\thighQualityHetCount\thighQualityLength\thighQualityHeterozygosity\tLOHlength\tannotationLohFilteredHetCount\tannotationLohFilteredLength\tannotationLohFilteredHeterozygosity\tdepthFilteredHetCount\tdepthFilteredLength\tdepthFilteredHeterozygosity\tsitesIn950kbHetCount\tsitesIn950kbLength\tsitesIn950kbHeterozygosity\n";

print TABLE "Strain & MTL & MLST Clade (FP\\footnote\\ ) & Heterozygosity\\footnote\\ & LOH length (Mbp)\\footnote\\ & Filtered length (Mbp)\\footnote\\ & Filtered heterozygosity\\footnote & Heterozy-gosity in 950 kb\\footnote\\ \\\\ \\hline\n";

my $sumoak950 = 0;
my $sumclinic950 = 0;

foreach my $het (reverse sort keys %hqh2strains) {
	my $strain = $hqh2strains{$het};

	my $mtl;
	if ($mtl{$strain} eq "het") { $mtl = "\$a/\\alpha\$"; }
	elsif ($mtl{$strain} eq "homa") { $mtl = "\$a/a\$"; }
	elsif ($mtl{$strain} eq "homalpha") { $mtl = "\$\\alpha/\\alpha\$"; }
	else { $mtl = $mtl{$strain}; } 

	if ($clades{$strain} eq "1") { $clades{$strain} = "1 (I)" }
	if ($clades{$strain} eq "2") { $clades{$strain} = "2 (II)" }
	if ($clades{$strain} eq "3") { $clades{$strain} = "3 (III)" }
	if ($clades{$strain} eq "4") { $clades{$strain} = "4 (SA)" }
	if ($clades{$strain} eq "6") { $clades{$strain} = "6 (I)" }
	if ($clades{$strain} eq "8") { $clades{$strain} = "8 (SA)" }
	if ($clades{$strain} eq "11") { $clades{$strain} = "11 (E)" }

	if ($strain =~/NCYC4/) { $sumoak950 += $allhps{$strain}/$allength; }
	elsif ($strain ne "1AA") { $sumclinic950 += $allhps{$strain}/$allength; }

	print STABLE "$strain\t$mtl{$strain}\t$clades{$strain}\t$hqps{$strain}\t$hql{$strain}\t$hqh{$strain}\t$lohl{$strain}\t$fps{$strain}\t$fl{$strain}\t$fh{$strain}\t$dfps{$strain}\t$dfl{$strain}\t$dfh{$strain}\t$allhps{$strain}\t$allength\t".($allhps{$strain}/$allength)."\n";

	print TABLE "$strain & $mtl & $clades{$strain} & ".sprintf("%.4f",$hqh{$strain})." & ".sprintf("%.1f",$lohl{$strain}/1000000)." & ".sprintf("%.1f",$dfl{$strain}/1000000)." & ".sprintf("%.4f",$dfh{$strain})." & ".sprintf("%.4f",$allhps{$strain}/$allength)." \\\\ ";

	if ($strain eq "NCYC4144") { 
		print TABLE "\n\\textbf{3 oak strains} &  & Mean: & ".sprintf("%.4f",$sumhqhncyc/$nncyc)." & ".sprintf("%.1f",($sumlohncyc/$nncyc)/1000000)." & &".sprintf("%.4f",$sumdfhncyc/$nncyc)." &  ".sprintf("%.4f",$sumoak950/$nncyc)."\\\\ \\\\ \n";
	}
	elsif ($strain eq "1AA") {	# is this the last entry?
		print TABLE "\n\\textbf{".($nstrains-4)." clinical strains}\\footnote\\ &  & Mean: & ".sprintf("%.4f",$sumhqhclinic/$nclinic)." & ".sprintf("%.1f",($sumlohclinic/$nclinic)/1000000)." & & ".sprintf("%.4f",$sumdfhclinic/$nclinic)." & ".sprintf("%.4f",$sumclinic950/$nclinic)." \\\\ \\hline \n";
	}
	else { print TABLE "\n"; }	

}

close TABLE;
close STABLE;


print "\nThere are summaries in $tableDH and $stable\n\n";

#################
## SUBROUTINES ##
#################

sub printheader {

	my $runname = shift;

	print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N $runname
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=168:00:00
#PBS -l mem=64gb

#PBS -M doudab\@uga.edu 
#PBS -m ae


cd $dir

module load R/3.2.3
";

}


sub printphyloheader {

	my $runname = shift;

	print SHOUT "

#PBS -S /bin/bash
#PBS -q batch
#PBS -N $runname
#PBS -l nodes=1:ppn=4:dbnode
#PBS -l walltime=168:00:00
#PBS -l mem=64gb

#PBS -M doudab\@uga.edu 
#PBS -m ae


cd $dir

module load raxml/8.1.20
module load R/3.2.3
";

	
}


# SUBROUTINE FOR PRINTING PLOTS FOR EACH STRAIN

sub printstrain {

	my $file = shift;
	my $strain;
	if ($file =~ /^(\S+?)\./m) { $strain = $1; }

	print RCMD "
data$strain<-read.table(\"$file\", header=T)
head(data$strain)
";
	
	my $seenchr = 0;
	foreach my $chr (@chr) {

	#	print "printing $file $chr\n";
		$seenchr++;

							# find out the length of each chromosome for adjusting plot widths with layout
		print RCMD "
# length of $strain $chr
max(data$strain\$pos[data$strain\$chr==\"$chr\"])
";
		if ($seenchr == 1) {			# the first entry for this strain
			print RCMD "
par(mar=c(2,6,1,0))
plot(c(1,max(data$strain\$pos[data$strain\$chr==\"$chr\"])),c(0,1),main=\"$chr\",type=\"n\",xlab=\'\',ylab=\'\')
title(ylab=\"Allele ratio\",line=2,cex.lab=1.5)
title(ylab=\"$strain\",line=4,cex.lab=1.5)
par(mar=c(2,0,1,0))
";
		}
		elsif ($seenchr == 8) {			# the last entry for this strain
			print RCMD "
par(mar=c(2,0,1,1))
plot(c(1,max(data$strain\$pos[data$strain\$chr==\"$chr\"])),c(0,1),main=\"$chr\",type=\"n\",yaxt=\'n\',xlab=\'\')
";			
		}
		else {
			print RCMD "
plot(c(1,max(data$strain\$pos[data$strain\$chr==\"$chr\"])),c(0,1),main=\"$chr\",type=\"n\",yaxt=\'n\',xlab=\'\')
";
		}
		annotate("$file$chr");

		print RCMD "
points(data$strain\$pos[data$strain\$chr==\"$chr\"],data$strain\$pALT[data$strain\$chr==\"$chr\"],main=\"$strain$chr\",pch=20,cex=0.1)
";
	}
}



# SUBROUTINE FOR PRINTING PLOTS FOR CHR OF EACH STRAIN

sub annotate {
						# annotate the plot with slightly transparent colored rectangles
	my $chr = shift;	# actually this is filechr

	my $height = 100000000;

	print RCMD "par(xpd=F)\n";	# do not print rectangles outside the plot	

	if (defined $astart{$chr}[0]) { 

		for (my $i=0; $i<@{$astart{$chr}}; $i++) {
			
#			print RCMD "rect($astart{$chr}[$i],0,$aend{$chr}[$i],$height,col=rainbow(1,alpha=0.3),border=rainbow(1,alpha=0.3))\n"; 

						# BITMAP TIFF DOES NOT SUPPORT TRANSPARENCY
			print RCMD "rect($astart{$chr}[$i],0,$aend{$chr}[$i],1,col=\"lightblue1\",border=\"lightblue1\")\n"; 		

		}
		print RCMD "par(xpd=T)\n";	# do print legend outside the plot	


#	for (my $j=1; $j<=@types; $j++) {
	#	print RCMD "legend(W[100],-0.15-($j/20),lty=1,legend=\"$types[$j-1]\",col=rainbow(".@types.",alpha=0.3)[$j],bty=\"n\")\n";	
#	}
		print RCMD "par(xpd=F)\n";	# don't print annotations outside the plot
	}
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

clade1col<-\"#dd1c77\"			# dark pink
clade6col<-\"#fa9fb5\"			# pink
clade4col<-\"#a7e0f2\" 			# light blue
clade8col<-\"cornsilk4\"
clade11col<-\"purple\"			
clade2col<-\"orange\"			
clade3col<-\"blue\"
clade5col<-\"red\"
clade7col<-\"#ffff99\"
clade9col<-\"#00441b\"
clade10col<-\"#a6cee3\"

wildcol<-\"#5aae61\"			# green
ncyccol<-\"brown\"
tc<-rep(\"brown\",length(t\$tip)) # default tip color = ncyccol = unknown


tc[grep(\"NCYC4\",t\$tip)]<-wildcol	
tc[grep(\"NCYC597\",t\$tip)]<-ncyccol	

tc[grep(\"1AA\",t\$tip)]<-clade1col
tc[grep(\"19F\",t\$tip)]<-clade1col
tc[grep(\"L26\",t\$tip)]<-clade1col
tc[grep(\"P37039\",t\$tip)]<-clade1col
tc[grep(\"12C\",t\$tip)]<-clade1col
tc[grep(\"P37005\",t\$tip)]<-clade1col
tc[grep(\"P37037\",t\$tip)]<-clade1col
tc[grep(\"P78048\",t\$tip)]<-clade1col
tc[grep(\"SC5314\",t\$tip)]<-clade1col


tc[grep(\"P94015\",t\$tip)]<-clade6col

tc[grep(\"GC75\",t\$tip)]<-clade4col
tc[grep(\"P75016\",t\$tip)]<-clade4col
tc[grep(\"P75063\",t\$tip)]<-clade4col
tc[grep(\"P87\",t\$tip)]<-clade4col


tc[grep(\"P60002\",t\$tip)]<-clade8col

tc[grep(\"P34048\",t\$tip)]<-clade3col
tc[grep(\"P78042\",t\$tip)]<-clade3col
tc[grep(\"P57055\",t\$tip)]<-clade3col

tc[grep(\"P57072\",t\$tip)]<-clade2col
tc[grep(\"P76055\",t\$tip)]<-clade2col
tc[grep(\"P76067\",t\$tip)]<-clade2col

tc[grep(\"P75010\",t\$tip)]<-clade11col


tc[grep(\"ST69\",t\$tip)]<-clade1col
tc[grep(\"ST79\",t\$tip)]<-clade1col
tc[grep(\"ST277\",t\$tip)]<-clade1col
tc[grep(\"ST103\",t\$tip)]<-clade1col
tc[grep(\"ST127\",t\$tip)]<-clade1col
tc[grep(\"ST408\",t\$tip)]<-clade1col
tc[grep(\"ST485\",t\$tip)]<-clade1col
tc[grep(\"ST1089\",t\$tip)]<-clade1col
tc[grep(\"ST1102\",t\$tip)]<-clade1col
tc[grep(\"ST1111\",t\$tip)]<-clade1col
tc[grep(\"ST1100\",t\$tip)]<-clade1col
tc[grep(\"ST1104\",t\$tip)]<-clade1col
tc[grep(\"ST1110\",t\$tip)]<-clade1col
tc[grep(\"ST1101\",t\$tip)]<-clade1col
tc[grep(\"ST347\",t\$tip)]<-clade5col
tc[grep(\"ST824\",t\$tip)]<-clade5col
tc[grep(\"ST886\",t\$tip)]<-clade5col
tc[grep(\"ST155\",t\$tip)]<-clade2col
tc[grep(\"ST232\",t\$tip)]<-clade2col
tc[grep(\"ST315\",t\$tip)]<-clade2col
tc[grep(\"ST538\",t\$tip)]<-clade11col
tc[grep(\"ST173\",t\$tip)]<-clade11col
tc[grep(\"ST514\",t\$tip)]<-clade11col
tc[grep(\"ST461\",t\$tip)]<-clade11col
tc[grep(\"ST1113\",t\$tip)]<-clade11col
tc[grep(\"ST186\",t\$tip)]<-clade3col
tc[grep(\"ST344\",t\$tip)]<-clade3col
tc[grep(\"ST726\",t\$tip)]<-clade3col
tc[grep(\"ST365\",t\$tip)]<-clade8col
tc[grep(\"ST90\",t\$tip)]<-clade8col
tc[grep(\"ST756\",t\$tip)]<-clade8col
tc[grep(\"ST729\",t\$tip)]<-clade8col
tc[grep(\"ST785\",t\$tip)]<-clade8col
tc[grep(\"ST1105\",t\$tip)]<-clade8col
tc[grep(\"ST1106\",t\$tip)]<-clade8col
tc[grep(\"ST1112\",t\$tip)]<-clade8col
tc[grep(\"ST409\",t\$tip)]<-clade6col
tc[grep(\"ST321\",t\$tip)]<-clade6col
tc[grep(\"ST377\",t\$tip)]<-clade6col
tc[grep(\"ST124\",t\$tip)]<-clade4col
tc[grep(\"ST95\",t\$tip)]<-clade4col
tc[grep(\"ST659\",t\$tip)]<-clade4col
tc[grep(\"ST1101\",t\$tip)]<-clade4col
tc[grep(\"ST840\",t\$tip)]<-clade7col
tc[grep(\"ST238\",t\$tip)]<-clade7col
tc[grep(\"ST651\",t\$tip)]<-clade7col
tc[grep(\"ST918\",t\$tip)]<-clade9col
tc[t\$tip==\"ST111\"]<-clade9col
tc[grep(\"ST735\",t\$tip)]<-clade9col
tc[grep(\"ST1103\",t\$tip)]<-clade9col
tc[grep(\"ST304\",t\$tip)]<-clade10col




#plot(t,tip.color=tc,cex=1.2,main=\"$chr\")
plot(t,tip.color=tc,cex=1.2)

bs<-as.numeric(t\$node)				# show only bootstraps >= 50%
bs[bs<50]<-NA
#nodelabels(bs,frame=\"none\",cex=0.5,bg=\"white\",adj = c(1.1,-0.2))
add.scale.bar(0.003,0)

#nodelabels(bs,frame=\"none\",bg=\"white\",adj = c(1.5,-0.3))
nodelabels(bs,frame=\"none\",bg=\"white\",adj = c(1.1,-0.2))

";

}


