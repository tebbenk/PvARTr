#!/usr/bin/perl -w                                                                                                                                             \
                                                                                                                                                                            
use strict;
use warnings;

open (IN, "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/sample_lists/redos.txt");
my @sample_names;
my $sorted;
my $unsorted_bam;

while(<IN>) {
#	print;
	chomp;
	push @sample_names, $_;
}

#print join " ", @sample_names;

for my $sample (@sample_names) {
	print "$sample\n";
	$unsorted_bam = "${sample}_pv_subset_rmdup.bam";
	print "$unsorted_bam\n";
	$sorted = "${sample}_pv_subset_rmdup_sorted.bam";
	
	system("/usr/local/packages/samtools-1.9/bin/samtools sort $unsorted_bam -o $sorted"); #sort bam file
}

exit;
