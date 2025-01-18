#!/usr/bin/perl -w                                                                                                                                             \
                                                                                                                                                                            
use strict;
use warnings;

open(MASK, "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/GATK_genomes/Pv4_mask.txt");
open(OUT, ">", "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/GATK_genomes/Pv4_mask_sequence.txt");

my $chr;
my $start;
my $end;
my $position;

while(my $line = <MASK>){
	my @info = split(/\t+/, $line);
	$chr = $info[0];
	$start = $info[1];
	$end = $info[2];
	for(my $x = $start; $x <= $end; $x++){
		$position = "$chr\_$x";
		print OUT "$position\n";
	}
}

exit;
