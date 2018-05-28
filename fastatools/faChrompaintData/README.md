# Usage Instructions for faChrompaint.pl with demo data files

With the following command, faChrompaint, will create plots for each alignment in files that begin with 'chr'. 
Each alignment will be broken into blocks of the default window size. 
Pairwise differences in each block are summarized in P34048.pDiffs.tsv and the nearest other strain is identified in P34048.nearest.tsv
each rectangle is colored according to the clade of the nearest sequence using the clade and color descriptions in 'clades' and 'colors'
The -M option means that the maximum difference allowed in order to give a clade assignment is 0.00066 (ie 0.066%) because 90% of within-clade differences for these C. albicans data are below this value

```
faChrompaint.pl -I chr -r P34048 -c clades -C colors -M 0.00066 -e 'NCYC4146 1AA SC5314_A'
```

To replot while only changing the colors of clades without the (slow) finding of most similar strains use the following command

```
faChrompaint.pl -p P34048.nearest.tsv -c clades -C colors 
```

NOTE: the -p plotting mode does not allow changes to which strains are excluded from the analysis or changes to any other parameters used when generating the .nearest.tsv file


