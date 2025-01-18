#!/usr/bin/perl -w                                                                                                                                             \
                                                                                                                                                                            
use strict;
use warnings;

my $genome = "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/genomes/PlasmoDB-54_PvivaxP01_Genome.fasta";
my $sorted_bam;

##Open list of sample names which are the directories and put into an array##
open (IN, "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/sample_lists/susceptible.txt");
my @sample_names;

while(<IN>) {
	chomp;
	push @sample_names, $_;
}

my $sample_directory;

for my $sample (@sample_names) {
	chomp $sample;
	print "$sample\n";
	my $bam = "$sample\_sorted.bam";
	my $bam_path = "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/WGS_ART/susceptible/sorted_bam_bowtie2/$bam";
#	print "$sorted_bam\n";
	my $mpileup = "$sample.mpileup";

	system("/usr/local/packages/samtools-1.9/bin/samtools mpileup -f $genome $bam_path > $mpileup"); 
}

close IN;

exit;