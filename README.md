# scripts

## Scripts used in Bensasson (2011) BMC Evolutionary Biology 11:211

### annotateCENs.pl

 A simple perl script to annotate Saccharomyces centromeres: This script uses a regular expression to recognise the consensus sequence motifs for CDEI and CDEIII described in Baker and Rogers 2005, Genetics 171(4):1463-1475.

### fourgamete.pl

A perl script implementing the four-gamete test for recombination described in Hudson and Kaplan (1985) Genetics 111(1)147-164. This test will not tell you whether or not there is recombination in your data. The script only outputs a summary of the alleles at each site, an alignment of only the segregating sites for easier visualisation of pairs of sites, and a list of the pairs of sites detected by the four-gamete test. In cases, where there are only 2 alleles at each site, these are cases where all 4 possible pairs of SNPs occur at the two sites. This implies either homoplasy (convergent mutation) or a recombination somewhere between the 2 sites. 

e.g. this combination can only be explained by homoplasy or recombination and fourgamete.pl would report positions 1 and 37 as a pair of sites detected by the four-gamete test
```
         position1 position37
taxon1   T         G
taxon2   T         A
taxon3   G         A
taxon4   G         G
```

For DNA sequence data, homoplasy is a likely explanation for many of these sites. This script does not summarise the data or tell you how many recombination events there are or where.

## Scripts used in Tilakaratna and Bensasson (unpublished)

### alcat.pl

A perl script that concatenates multiple alignment files in fasta format into a single large alignment file

### structureInfile.pl

Converts DNA sequence alignments in fasta format to STRUCTURE input files that summarises bases at variable sites.

### structureShell.pl

Runs STRUCTURE one time for each value of K in a specified range (e.g. from 1 to 10).

### structurePrint.pl

Plots STRUCTURE results as barplots using R with user control of colors.


