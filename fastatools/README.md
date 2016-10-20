## fastatools/

This directory contains scripts used to manipulate DNA sequence in fasta format:

### alcat.pl

Used in Tilakaratna and Bensasson (unpublished).

A perl script that concatenates multiple alignment files in fasta format into a single large alignment file

### fa2phylip.pl

This perl script converts fasta format sequence into alignments in phylip format.

### faChoose.pl

A perl script to choose a subset of sequences from a fasta file. Currently, the script searches for the names provided by the user in a way that is case insensitive, and the pattern can be matched anywhere in the first word of the fasta name descriptor line.

### fastaLength.pl

A perl script to summarise the length of DNA sequences in a fasta file. The -g option is useful for showing the ungapped length of DNA sequences in an alignment.


