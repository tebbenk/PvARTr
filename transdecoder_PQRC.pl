#!/usr/bin/perl -w                                                                                                                                             \
                                                                                                                                                                            
use strict;
use warnings;

my @fastafiles = </local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/WGS_ART/all_vcf/ubp1redo/*.fa>;
my $transcripts;

for my $fasta (@fastafiles) {
	my @files = split /\//, $fasta;
	$transcripts = $files[9];	
	#print "$transcripts\n";
	
	system("/usr/local/packages/transdecoder-5.7.1/TransDecoder.LongOrfs -t $transcripts");
	
	system("/usr/local/packages/transdecoder-5.7.1/TransDecoder.Predict -t $transcripts --no_refine_starts");

}

exit;
