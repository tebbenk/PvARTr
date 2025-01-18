#!/usr/bin/perl -w                                                                                                                                             \

use strict;
use warnings;

my @bamfiles = </local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_longitudinal/bams_71924/bam_human/*RNA.bam>;
my $bam;
my $bam_assigned;

open(OUT, ">", "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_longitudinal/bams_71924/bam_human/mapped_human_82124.txt");
print OUT "Sample\tTotal\tMappedHumanReads\tPercent\n";

foreach $bam (@bamfiles){
open (BAM,"/usr/local/packages/samtools-1.9/bin/samtools view -h $bam |");
my @files = split /\//, $bam;
my $sample_name = $files[8];	
(my $sample = $sample_name) =~ s/\.[^.]+$//;
print "$sample\n";
my $bam_human = "${sample}_subset_human.bam";
my $total =0;
my $human = 0;
my $percent = 0;

open(HUMAN, "| /usr/local/packages/samtools-1.9/bin/samtools view -Sb - > $bam_human");

	while( my $line = <BAM>){
    		chomp $line;
    		if($line =~ m/^@/) { 
    			print HUMAN "$line\n";
    			next;
	  		}
    	my @sam = split(/\t+/, $line);  ## splitting SAM line into array
    	if($sam[1] == 83 or $sam[1] == 163 or $sam[1] == 99 or $sam[1] == 147){
	   		$total++;
			if($sam[2] =~ "N[C,T,W]_"){
				print HUMAN "$line\n";
     			$human++;
     		}
     	}    	 
	}
	$percent = ($human/$total) * 100;
	print OUT "$sample";
	print OUT "\t$total";
	print OUT "\t$human";
	print OUT "$percent";
	print OUT "\n";
	close HUMAN;
}

close OUT; 

exit; 