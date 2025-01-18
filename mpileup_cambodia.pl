#!/usr/bin/perl -w                                                                                                                                             \
                                                                                                                                                                            
use strict;
use warnings;

open (IN, "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/sample_lists/redos.txt");
my @sample_names;
my $sorted_bam;
my $mpileup;

while(<IN>) {
	chomp;
	push @sample_names, $_;
}

for my $sample (@sample_names) {
	print "$sample\n";
	$sorted_bam = "${sample}_pv_subset_rmdup_sorted.bam";
	$mpileup = "${sample}.mpileup";

#	print "$sorted_bam\n";
#	print "$mpileup\n";
	my $genome = "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/genomes/PlasmoDB-54_PvivaxP01_Genome.fasta";
	
	system("/usr/local/packages/samtools-1.9/bin/samtools mpileup -f $genome $sorted_bam > $mpileup");
	
}

exit;
