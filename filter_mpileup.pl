#!/usr/bin/perl -w                                                                                                                                             \
                                                                                                                                                                            
use strict;
use warnings;

my %keep;

#Read in file containing all locations in the genome EXCEPT multigene families - you should have access to this file
open(GENES, "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/genomes/pv_genes_multigenefams_removed.txt");

#Puts all regions you do want (excludes multi-gene families) in a hash
print "generating hash of genes to keep\n";
while(my $line = <GENES>){
#	print $line;
	my @info = split(/\t+/, $line);
#	print "$info[1]\n";
	my $chr = $info[0];
#	print "$chr\n";
	my $start = $info[1];
	my $end = $info[2];
#	print "$chr\t$start\t$end\n";
	for(my $i = $start; $i <= $end; $i++){
		my $position = $chr."_$i";
		$keep{"$position"} = 1;
#		print "$position\n";
	}
}

close(GENES);

#Read in your mpileup files - you just need to change the path to whatever folder your mpileup files are in. *.mpileup will read in any file into the array that ends in .mpileup
my @mpileup_files = </local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/WGS_ART/resistant/sorted_bams_bowtie2/mpileup/*.mpileup>;
my $mpileup;

#Loops through each mpileup file that you have and 
for $mpileup (@mpileup_files){
	print "$mpileup\n";
	
	my @files = split /\//, $mpileup;
	
	#Change the number to whatever position in the path corresponds to ".mpileup"; it might be different based on how many sub folders you have
	my $sample_name = $files[10];	
	(my $sample = $sample_name) =~ s/\.mpileup//;
	print "$sample\n";
	
	#Generates a new file that ends in filtered.mpileup for each input mpileup
	my $out = $sample."_filtered.mpileup";
	
	#opens that file for writing - change the path to whatever folder you want to write the files to
	open(OUT, ">", "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/WGS_ART/resistant/sorted_bams_bowtie2/mpileup/filtered_mpileup/$out");
	
	#opens each mpileup
	open(IN, $mpileup);
	
	#reads through each line of the mpileup and prints only lines that do not correspond to multigene families	
	while(my $pileup_line = <IN>){
#		print "$pileup_line\n";
		my @info = split(/\t+/, $pileup_line);
		my $chr = $info[0];
		my $loc = $info[1];
		my $position = $chr."_$loc";
#		print "$position\n";
		if (defined $keep{$position}){
#			print "$position\n";
#			print "$pileup_line\n";
			print OUT "$pileup_line";
		}
	}
	close(OUT);
}


exit;