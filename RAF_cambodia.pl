#!/usr/bin/perl -w                                                                                                                                             \
                                                                                                                                                                            
use strict;
use warnings;

open (SAMPLES, "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/sample_lists/sample_69.txt");
my @sample_names;

while(<SAMPLES>) {
	chomp;
	push @sample_names, $_;
}

my $mpileup;
my @read;
my $position;
my $allele;
my $bases;
my $reference;
my $RAF; 

for my $sample (@sample_names){
	print "$sample\n";
	my $out = $sample."_RAF.txt";
	
	open(OUT, ">", "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/WGS_ART/resistant/sorted_bams_bowtie2/mpileup/filtered_mpileup/RAF/$out");
	print OUT "chromosome\tposition\tcoverage\treads_with_ref_allele\tRAF\n";
	
	chomp $sample;
	my $mpileup_path = "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/WGS_ART/resistant/sorted_bams_bowtie2/mpileup/filtered_mpileup/${sample}_filtered.mpileup";
	open(IN, $mpileup_path);
	while(my $line = <IN>){
		chomp $line;
#		print $line;
		$reference = 0;
		$bases = 0;
		@read = split(/\t+/, $line);  ## splitting line into array
		next if $read[4] =~ m/\+/;
		for($position = 0; $position < length $read[4]; $position++){
			$allele = substr($read[4], $position, 1);
			if($allele eq "\." or $allele eq "\," or $allele eq "\>" or $allele eq "\<" or $allele eq "A" or $allele eq "a" or $allele eq "C" or $allele eq "c" or $allele eq "G" or $allele eq "g" or $allele eq "T" or $allele eq "t"){
				$bases++;				
			}
			if($allele eq "\." or $allele eq "\," or $allele eq "\>" or $allele eq "\<"){
				$reference++;
			}	
		}
		if($bases >= 50){
			$RAF = $reference/$bases;
 			print OUT "$read[0]\t";
 			print OUT "$read[1]\t";
 			print OUT "$read[3]\t";
 			print OUT "$reference\t";
 			print OUT "$RAF\n";
		}	
	}
	
	close(OUT);
}

exit;