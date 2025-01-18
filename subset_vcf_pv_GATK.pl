#!/usr/bin/perl -w                                                                                                                                             \

##Filters our variable genome regions, removes indels

                                                                                                                                                                            
use strict;
use warnings;

my %mask;

open(GENES, "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/GATK_genomes/Pv4_mask.txt");

while(my $line = <GENES>){
	my @info = split(/\t+/, $line);
	my $chr = $info[0];
	my $start = $info[1];
	my $end = $info[2];
	for(my $i = $start; $i <= $end; $i++){
		my $position = $chr."_$i";
		$mask{"$position"} = 1;
	}
}

close(GENES);

my $vcf = "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/WGS_ART/resistant/GATK_files/resistant_51424_filt2.vcf.recode.vcf";

my $out = "resistant_moimix_ready_51524.vcf";
#	print "$out\n";
open(OUT, ">", "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/WGS_ART/resistant/GATK_files/$out");

open(IN, $vcf);
	
while(my $vcf_line = <IN>){
	chomp $vcf_line;
	if($vcf_line =~ m/^#/) {
		if($vcf_line =~ "N[CTW]_"){
			next;
		} 
		else{
			print OUT "$vcf_line\n";
			next;
		}
	}
	my @info = split(/\t+/, $vcf_line);
	my $type = $info[7];
	my $chr = $info[0];
 	my $loc = $info[1];
 	my $position = $chr."_$loc";
 	if (defined $mask{$position}){
		next;
 	}
 	else{
 	 	if ($type =~ "INDEL"){
 			next;
 		}
 		else{
 			print OUT "$vcf_line\n";
 		}
 	}
}
close(OUT);



exit;