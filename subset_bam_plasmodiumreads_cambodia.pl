#!/usr/bin/perl -w                                                                                                                                             \

use strict;
use warnings;

my @bamfiles = </local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_longitudinal/*RNA.bam>;
my $bam;
my $bam_assigned;

open(OUT, ">", "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_longitudinal/mappingstats_pv_005_D0.txt");
print OUT "Sample\tTotal\tMappedPvReads\tPercentPv\n";

foreach $bam (@bamfiles){
	open (BAM,"/usr/local/packages/samtools-1.9/bin/samtools view -h $bam |");
	my @files = split /\//, $bam;
	my $sample_name = $files[6];	
	(my $sample = $sample_name) =~ s/\.[^.]+$//;
	print "$sample\n";
	my $bam_pv = "${sample}_subset_pv.bam";
	my $total = 1;
	my $pv = 0;
	my $percent;

	open(PV, "| /usr/local/packages/samtools-1.9/bin/samtools view -Sb - > $bam_pv");	

	while( my $line = <BAM>){
    	chomp $line;
    	if($line =~ m/^@/) { 
			print PV "$line\n";
  			next;
	  		}
    	my @sam = split(/\t+/, $line);  ## splitting SAM line into array
    	if($sam[1] == 83 or $sam[1] == 163 or $sam[1] == 99 or $sam[1] == 147){
	   		$total++;
      		if($sam[2] =~ "Pv"){
      			print PV "$line\n";
      			$pv++;
 	     	}
     	}    	 
	}
	$percent = ($pv/$total) * 100;
	print OUT $sample;
	print OUT "\t$total";
	print OUT "\t$pv";
	print OUT "\t$percent";
	print OUT "\n";
	close PV;
	close BAM;
}

close OUT; 

exit; 