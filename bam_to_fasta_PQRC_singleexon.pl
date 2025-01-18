#!/usr/bin/perl -w                                                                                                                                             \
                                                                                                                                                                            
use strict;
use warnings;

my @bamfiles = </local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/WGS_ART/susceptible/sorted_bams_bowtie2/*.bam>;
my $genome = "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/genomes/PlasmoDB-54_PvivaxP01_Genome.fasta";
my $K13;
my $MRP2;
my $PIP3;
my $AP2Mu;
my $VP3;
my $Falcilysin;

for my $bam (@bamfiles) {
	my @files = split /\//, $bam;
	my $input_bam = $files[9];
	(my $sample = $input_bam) =~ s/_sorted.bam//;
	print "$sample\n";
	
	
	my $vcf = "${sample}.vcf.gz";
	my $full_fasta = "${sample}.fa";
	$K13 = "${sample}_K13.fa";
	$MRP2 = "${sample}_MRP2.fa";
	$PIP3 = "${sample}_PIP3.fa";
	$AP2Mu = "${sample}_AP2Mu.fa";
	$VP3 = "${sample}_VP3.fa";
	$Falcilysin = "${sample}_falcilysin.fa";
	
	
	system("bcftools mpileup -Ou -f $genome $input_bam | bcftools call -Ou -mv --ploidy 1 | bcftools norm -f $genome -Oz -o $vcf");
	
	system("tabix $vcf");
	
	system("bcftools consensus -f $genome $vcf > $full_fasta");
	
	#K13
	system("samtools faidx $genome PvP01_12_v2:482687-487553 | bcftools consensus $vcf -o $K13");
	
	#MRP2
	system("samtools faidx $genome PvP01_14_v2:2051166-2059988 | bcftools consensus $vcf -o $MRP2");
	
	#PIP3
	system("samtools faidx $genome PvP01_10_v2:826183-831486 | bcftools consensus $vcf -o $PIP3");
	
	#AP2Mu
	system("samtools faidx $genome PvP01_14_v2:1556502-1561611 | bcftools consensus $vcf -o $AP2Mu");
	
	#VP3
	system("samtools faidx $genome PvP01_09_v2:725627-729307 | bcftools consensus $vcf -o $VP3");
	
	#Faliclysin
	system("samtools faidx $genome PvP01_11_v2:451243-458605 | bcftools consensus $vcf -o $Falcilysin");
	
}

exit;