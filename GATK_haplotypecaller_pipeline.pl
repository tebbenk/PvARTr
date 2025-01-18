#!/usr/bin/perl -w                                                                                                                                             \
                                                                                                                                                                            
use strict;
use warnings;


#Mapping through GATK HaplotypeCaller to generate single sample VCF files
#Important that $ref is indexed according to: https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format
#Important that $known is indexed according to: https://gatk.broadinstitute.org/hc/en-us/articles/360037594851-IndexFeatureFile
#Must load python and java before running


#my $genome = "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/GATK_genomes/pv_bowtie_index_GATK";
#my $known = "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/GATK_genomes/Pv4_mask_regions.bed";
my $ref = "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/GATK_genomes/PvivaxP01.genome.fasta";
my $sam;
my $log;
my $bam; 
my $sorted;
my $dedup;
my $dedup_met;
my $table;
my $recal;
my $gvcf;

##Open list of sample names which are the directories and put into an array##
open (IN, "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/sample_lists/resistant.txt");
my @sample_names;

while(<IN>) {
	chomp;
	push @sample_names, $_;
}

my $sample_directory;

for my $sample (@sample_names) {
	chomp $sample;
	print "$sample\n";
#	$sample_directory = "/local/projects-t4/aberdeen2ro/SerreDLab-4/raw_reads/2024-03-05_UMB/$sample/ILLUMINA_DATA";
#	print "$sample_directory\n";
#	sleep(1);
#	$sam = "${sample}.sam";
#	$log = "${sample}_log.txt";
#	$bam = "${sample}_bowtie.bam";
	$sorted = "${sample}_bowtie_sorted.bam";
#	$dedup = "${sample}_bowtie_sorted_dedup.bam";
#	$dedup_met = "${sample}_dedup_metrics.txt";
#	$table = "${sample}_recal.table";
#	$recal = "${sample}_recal.bam";
	$gvcf = "${sample}.g.vcf";
#	print "$sam\n$log\n";
#	sleep(1);

# #Open directory for each sample and pull out READ 1 files ##		
# 	opendir(DIR, $sample_directory);
# 	my @R1files = grep(/R1.fastq.gz$/,readdir(DIR));
# 	closedir(DIR); 
# 		
# 	my $R1file;
# 	my $R1_filepath;
# 		
# 	foreach $R1file (@R1files) {
# 		$R1_filepath = "$sample_directory/${R1file}";
# # 		print "$R1_filepath\n";
# 	}
# 
# ##Open directory for each sample and pull out READ 2 files ##
# 	opendir(DIR, $sample_directory);
# 	my @R2files = grep(/R2.fastq.gz$/,readdir(DIR));
# 	closedir(DIR);
# 		
# 	my $R2file;
# 	my $R2_filepath;
# 		
# 	foreach $R2file (@R2files) {
# 		$R2_filepath = "$sample_directory/${R2file}";
# #		print "$R2_filepath\n";
# 	}
	
	#Map
# 	print "mapping\n";
# 	system("/usr/local/packages/bowtie2-2.4.4/bowtie2 -x $genome -1 $R1_filepath -2 $R2_filepath -S $sam --rg-id $sample --rg SM:$sample --rg PL:Illumina"); 
# 	
# 	#Convert sam to bam
# 	print "converting to bam\n";
# 	system("/usr/local/packages/samtools-1.9/bin/samtools view -S -b -F 4 $sam > $bam"); 
# 
# 	#Remove sam file
# 	print "deleting sam\n";
# 	system("unlink $sam");
# 	
# 	#Sort bam using GATK sort function 
# 	print "sorting\n";
# 	system("/usr/local/packages/gatk-4.2.2.0/gatk SortSam -I $bam -O $sorted -VALIDATION_STRINGENCY LENIENT -USE_JDK_DEFLATER true -USE_JDK_INFLATER true --SORT_ORDER coordinate");
# 	
# 	#Mark duplicates
# 	print "Marking duplicates\n";
# 	system("/usr/local/packages/gatk-4.2.2.0/gatk MarkDuplicates --REMOVE_DUPLICATES --USE_JDK_DEFLATER true --USE_JDK_INFLATER true -I $sorted -O $dedup --METRICS_FILE $dedup_met");
# 	
# 	#Base recalibration
# 	print "Base Recalibration\n";
# 	system("/usr/local/packages/gatk-4.2.2.0/gatk BaseRecalibrator -I $dedup -O $table -R $ref --known-sites $known");
# 	
# 	#Apply base recalibration
# 	print "Applying base recalibration\n";
# 	system("/usr/local/packages/gatk-4.2.2.0/gatk ApplyBQSR -I $dedup -O $recal -R $ref --bqsr-recal-file $table");
	
	#Index
#	print "Indexing\n";
#	system("/usr/local/packages/samtools-1.9/bin/samtools index $sorted");
	
	#Haplotypecaller
	print "Running HaplotypeCaller";
	system("/usr/local/packages/gatk-4.2.2.0/gatk HaplotypeCaller -I $sorted -O $gvcf -R $ref -ERC GVCF");

}

close IN;