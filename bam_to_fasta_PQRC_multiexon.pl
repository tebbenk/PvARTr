#!/usr/bin/perl -w                                                                                                                                             \
                                                                                                                                                                            
use strict;
use warnings;

my @vcf_files = </local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/WGS_ART/all_vcf/*.vcf.gz>;
my $genome = "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/genomes/PlasmoDB-54_PvivaxP01_Genome.fasta";
# my $EPS15_1;
# my $EPS15_2;
# my $EPS15_3;
# my $EPS15_4;
# my $EPS15_5;
# my $EPS15_6;
# my $EPS15_7;
# my $EPS15_complete;
# my $EPS15_complete_merged;
my $UBP1_1;
my $UBP1_2;
my $UBP1_3;
my $UBP1_complete;
my $UBP1_complete_merged;
# my $coronin_1;
# my $coronin_2;
# my $coronin_3;
# my $coronin_complete;
# my $coronin_complete_merged;
# my $PATPL1_1;
# my $PATPL1_2;
# my $PATPL1_3;
# my $PATPL1_complete;
# my $PATPL1_complete_merged;


for my $vcf (@vcf_files) {
	my @files = split /\//, $vcf;
#	my $vcf = $files[10];
	(my $sample = $files[8]) =~ s/.vcf.gz//;
#	print "$vcf\n";
	
	
#	my $vcf = "${sample}.vcf.gz";
# 	$EPS15_1 = "${sample}_EPS15_exon1.fa";
# 	$EPS15_2 = "${sample}_EPS15_exon2.fa";
# 	$EPS15_3 = "${sample}_EPS15_exon3.fa";
# 	$EPS15_4 = "${sample}_EPS15_exon4.fa";
# 	$EPS15_5 = "${sample}_EPS15_exon5.fa";
# 	$EPS15_6 = "${sample}_EPS15_exon6.fa";
# 	$EPS15_7 = "${sample}_EPS15_exon7.fa";
# 	$EPS15_complete = "${sample}_EPS15_concat.fa";
# 	$EPS15_complete_merged = "${sample}_EPS15_complete.fa";
	
	$UBP1_1 = "${sample}_UBP1_exon1.fa";
	$UBP1_2 = "${sample}_UBP1_exon2.fa";
	$UBP1_3 = "${sample}_UBP1_exon3.fa";
	$UBP1_complete = "${sample}_UBP1_concat.fa";
	$UBP1_complete_merged = "${sample}_UBP1_complete.fa";
	
# 	$coronin_1 = "${sample}_coronin_exon1.fa";
# 	$coronin_2 = "${sample}_coronin_exon2.fa";
# 	$coronin_3 = "${sample}_coronin_exon3.fa";
# 	$coronin_complete = "${sample}_coronin_concat.fa";
# 	$coronin_complete_merged = "${sample}_coronin_complete.fa";
# 	
# 	$PATPL1_1 = "${sample}_PATPL1_exon1.fa";
# 	$PATPL1_2 = "${sample}_PATPL1_exon2.fa";
# 	$PATPL1_3 = "${sample}_PATPL1_exon3.fa";
# 	$PATPL1_complete = "${sample}_PATPL1_concat.fa";
# 	$PATPL1_complete_merged = "${sample}_PATPL1_complete.fa";
# 
# 
# 	#EPS15
# 	system("samtools faidx $genome PvP01_06_v2:447172-447787 | bcftools consensus $vcf -o $EPS15_1");
# 	system("samtools faidx $genome PvP01_06_v2:448015-448191 | bcftools consensus $vcf -o $EPS15_2");
# 	system("samtools faidx $genome PvP01_06_v2:448343-450195 | bcftools consensus $vcf -o $EPS15_3");
# 	system("samtools faidx $genome PvP01_06_v2:450367-450759 | bcftools consensus $vcf -o $EPS15_4");
# 	system("samtools faidx $genome PvP01_06_v2:450942-451065 | bcftools consensus $vcf -o $EPS15_5");
# 	system("samtools faidx $genome PvP01_06_v2:451221-451314 | bcftools consensus $vcf -o $EPS15_6");
# 	system("samtools faidx $genome PvP01_06_v2:451487-453922 | bcftools consensus $vcf -o $EPS15_7");
# 	system("awk 'NR==1 || FNR > 1' *EPS15_exon* > $EPS15_complete");
# 	system("seqkit -w 0 seq $EPS15_complete -o $EPS15_complete_merged");
# 	
	#UBP1_1
	system("samtools faidx $genome PvP01_02_v2:413850-414192 | bcftools consensus $vcf -o $UBP1_1");
	system("samtools faidx $genome PvP01_02_v2:414398-417125 | bcftools consensus $vcf -o $UBP1_2");
	system("samtools faidx $genome PvP01_02_v2:417242-425258 | bcftools consensus $vcf -o $UBP1_3");
	system("awk 'NR==1 || FNR > 1' *UBP1_exon* > $UBP1_complete");
	system("seqkit -w 0 seq $UBP1_complete -o $UBP1_complete_merged");
	
# 	#Coronin
# 	system("samtools faidx $genome PvP01_14_v2:2907698-2907731 | bcftools consensus $vcf -o $coronin_1");
# 	system("samtools faidx $genome PvP01_14_v2:2907871-2907962 | bcftools consensus $vcf -o $coronin_2");
# 	system("samtools faidx $genome PvP01_14_v2:2908226-2909797 | bcftools consensus $vcf -o $coronin_3");
# 	system("awk 'NR==1 || FNR > 1' *coronin_exon* > $coronin_complete");
# 	system("seqkit -w 0 seq $coronin_complete -o $coronin_complete_merged");
# 	
# 	#PATPL1
# 	system("samtools faidx $genome PvP01_04_v2:631410-632744 | bcftools consensus $vcf -o $PATPL1_1");
# 	system("samtools faidx $genome PvP01_04_v2:632965-633725 | bcftools consensus $vcf -o $PATPL1_2");
# 	system("samtools faidx $genome PvP01_04_v2:633878-635681 | bcftools consensus $vcf -o $PATPL1_3");
# 	system("awk 'NR==1 || FNR > 1' *PATPL1_exon* > $PATPL1_complete");
# 	system("seqkit -w 0 seq $PATPL1_complete -o $PATPL1_complete_merged");
# 	
# 	
	#Remove all intermediate files
	system("rm *exon*");
	system("rm *concat.fa");
	

	
}

exit;