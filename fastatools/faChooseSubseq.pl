#!/usr/bin/perl

$program = 'faChooseSubseq.pl';
$format = 'unknown';
($outfile = $program) =~ s/(.*?).pl/$1.out/; 
($missfile = $program) =~ s/(.*?)\.pl/$1.missed/; 

use Getopt::Std;
getopts('i:n:o:m:s:pr',\%parameters);


if (exists($parameters{"i"})) { $seqfile = $parameters{"i"}; }
if (exists($parameters{"n"})) { $namefile = $parameters{"n"}; }
if (exists($parameters{"o"})) { $outfile = $parameters{"o"}; }
if (exists($parameters{"m"})) { $missfile = $parameters{"m"}; }


unless (exists $parameters{"i"}) {
	print "\nUSAGE: $program -i <FASTAfile> -n <names or NAMEfile> \n\n"; 
	print "  -n\tnames: command line in quotes, in blast output or list file\n";
	print "  -o\tname for file of chosen seqs\n";
	print "  -m\tname for file for missing names\n";
	print "  -p\tfor pattern match (name is not exact, slower?)\n";
	print "  -s\tsubstring coordinates e.g. '30..200'\n";
	print "  -r\trevcomp/-ve strand (substring will be of -ve)\n\n";

	exit;
}

print "\n\t< $program > - Getting seqs from a FASTAfile\n\n";



open NAMES, "<$namefile" or $names = $namefile; 
if (defined $names) {
	
	while ($names =~ /(\S+)/g) { push (@names, $1); } 	# recognising > 1 name within quotes ''
	
	foreach $name (@names) { 
		$wantedname{$name}++;
#		if ($name =~ /^(\w+?\.[a-z]\S+)/mi) { $format = 'read'; }
		if ($name =~ /^(\d+)$/m) { $uid_length = length($1); $format = 'UIDnumber'; }
		else { $format = 'unknown'; print "NB did not recognise the name format of $name\n"; }
	}
	print "will get $format @names from $seqfile\n";
}
else {
	print "\nreading $namefile.\n";
#	while (<NAMES>) {
#		if (/^(\w+?\.[a-z]\S+?)/mi) { $format = 'read';
#		}
#	}
	unless ($format eq 'reads') {
		open NAMES, "<$namefile" or $names = $namefile; 
		while (<NAMES>) {
			if (/^(\d+)$/mi) { 
				$uid_length = length($1);
				$format = 'UIDnumber';
				unless ($wantedname{$1}++) { push(@names, $1) }; 
			}
		}
		unless ($format eq 'UIDnumber') {
			open NAMES, "<$namefile" or $names = $namefile; 
			while (<NAMES>) {
				if (/^(\S+)/mi) { 
					$uid_length = length($1);
					unless ($wantedname{$1}++) { push(@names, $1) }; 
				}
			}
		}
	}
	print "Found ".@names." names in $format format\n";
}



open FASTA, "<$seqfile" or die "couldn't open $seqfile: $!";

print "\nreading $seqfile. Progress in millions of lines:\n\t"; 	
my $progress = 0; 

if ($format eq 'read') {
	while (<FASTA>) {
		$progress++;
		if ($progress > 1000000) { print "*"; $progress = 0; }

		if (/^>(\w+?\.[a-z]\S+?)(\s*.*?)$/mi) {
			$name = $1; 
			$line{$name} = $1.$2;
			$line{$name} = chomp($line{$name});
			#	print "check: [$2]\n";		

			if (exists $wantedname{$1}) { push(@foundname, $name); }

		}

		elsif ((exists $wantedname{$name}) && (/^([acgtnx\?-]+)$/mi)) { $seq{$name} .= "$1"; }	# FASTA DNA SEQ
		elsif ((exists $wantedname{$name}) && (/^(\d[\d ]+)$/m)) { $seq{$name} .= "$1"; }	# FASTA QUAL
	}
#	die;
}

elsif (($format eq 'UIDnumber') || ($format eq 'unknown')) {
	while (<FASTA>) {
		$progress++;
		if ($progress > 1000000) { print "*"; $progress = 0; }
	
		if (/^>(\S+)(.*?)$/m) {
			$name = $1;
			$line{$name} = $1.$2;
			
			if (@foundname == @names) { last; }
			elsif (exists $wantedname{$name}) { 
				print "Found $name\n"; 
				push(@foundname, $name); 
				$found{$name}++;
			} 
			elsif (exists $parameters{"p"}) {
				foreach $wanted (@names) {
					if ($line{$name} =~ /$wanted/i) {
						print "Found $line{$name} : $name\n"; 
						push(@foundname, $name); 
						$found{$name}++;
					}
				}
			}
		} 
		elsif ((exists $found{$name}) && (/^([0-9a-z\?-]+)(.*?)$/mi)) { 
			$seq{$name} .= "$1"; 								# FASTA DNA SEQ
			my $extra = $2;			
			if (length($extra) > 0) { print "\nWARNING. ".length($extra)." character(s) omitted from sequence:\n[$extra]\n"; }
			if (/^([b,d,e-f,h-j]+)$/m) { $format = "staden_qual"; }				# FASTA STADEN QUAL

		}
		elsif ((exists $found{$name}) && (/^([A-Z-]+)$/m)) { $seq{$name} .= "$1"; }		# FASTA PROT SEQ
		elsif ((exists $found{$name}) && (/^(\d[\d ]+)$/m)) { $seq{$name} .= "$1"; $format = "qual"; }	# FASTA QUAL
		elsif (exists $found{$name}) { print "unrecognised seq: [$_]\n"; }
	}
}

print "Found ".@foundname." sequences. ";
if (@foundname > 0) { 
	if (exists $parameters{'a'}) { open OUT, ">>$outfile" or die "couldn't open $outfile: $!"; }
	else { open OUT, ">$outfile" or die "couldn't open $outfile: $!"; }
	foreach $name (@foundname) { 
		my $printname = $line{$name};
		
		
		if (exists $parameters{"s"}) {						# SUBSEQ (incl. >1 region)
			if ($parameters{"s"} =~ /\d+\.\.\d+/) {				# NOTE: REVCOMP AFTER SUBSTRING
				my $subseq;
				print "\n\nlength of $name : ".length($seq{$name})."bp\n";
				while ($parameters{"s"} =~ /(-?)(\d+)\.\.(\d+)/g) {
					my $start = $2;
					my $end = $3;
					if ($1 eq '-') { $start = 1; } 
					if (exists $parameters{"r"}) { $printname .= ":$start-$end"; }
					else { $printname = $name.":$start-$end"; }
					unless ($start == 0) {
						$subseq .= substr($seq{$name}, ($start - 1), ($end - $start + 1));
					}
					else {	$subseq .= substr($seq{$name}, ($start), ($end - $start)); }
					print "\t$start .. $end\t";
				}
				$seq{$name} = $subseq;
				print "printed : ".length($seq{$name})."bp\n\n";
				
			}
			else { print "Do not understand the subseq notation ($parameters{'s'})\n"; }
		}

		if (exists $parameters{"r"}) {						# REVCOMP (-ve strand)
			my $rev = reverse $seq{$name};					# NOTE: REVCOMP AFTER SUBSTRING
			($seq{$name} = $rev) =~ tr/acgtACGTmrwsykvhdbMRWSYKVHDB/tgcaTGCAkywsrmbdhvKYWSRMBDHV/;
			$printname .= "(-)";
		}

		print OUT ">$printname\n$seq{$name}\n"; 
	} 
	close OUT;
	print "Printed them to $outfile";
}

foreach $name (@names) {
	unless (exists $seq{$name}) { unless (exists $parameters{"p"}) {  push (@misses, $name); } }
}

print "\n".@misses." sequences were missing:\n";
if (@misses > 0) { 
	open MISSES, ">$missfile" or die "couldn't open $missfile: $!";
	foreach $miss (@misses) { print MISSES "$miss\n"; } 
	close MISSES;
	print "Printed them to $missfile.\n";
}
