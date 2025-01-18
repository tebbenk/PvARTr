#!/usr/bin/perl -w                                                                                                                                             \
                                                                                                                                                                            
use strict;
use warnings;

my %var_pos;
my $genotype;
my $allele;
my @read;
my $RAF; 
my $position;


#Generate hash of positions to keep
open(VCF, "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/WGS_ART/all_vcf/merged_allsamples_noMGF_filtered_onlychr_noapi_noMIT.vcf");
print "Generating hash\n";
while(my $line = <VCF>){	
	if($line =~ "#"){
		next;
	}
	else{
		my @info = split(/\t+/, $line);
		my $chr = $info[0];
		my $position = $info[1];
		my $chr_pos = $chr."_$position";
		$var_pos{"$chr_pos"} = 1;
	}	
}

my @mpileup_files = </local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/WGS_ART/all_vcf/snps_fastas/*snps.mpileup>;
open(OUT, ">", "/local/projects-t3/SerreDLab-3/kieran.tebben/PQRC_bulk/WGS_ART/all_vcf/snps_fastas/allsamples_variable_positions.fasta");

my %mpileup_positions;

for my $mpileup (@mpileup_files){
	open(IN, $mpileup);
	my @files = split /\//, $mpileup;
	my $sample_name = $files[9];	
	(my $sample = $sample_name) =~ s/_snps.mpileup//;
	print "$sample\n";
	print OUT ">$sample\n";
	
	#Put all positions into a hash
	while(my $pileup_line = <IN>){
		chomp $pileup_line;
		my $reference = 0;
		my $A = 0;
		my $T = 0;
		my $C = 0;
		my $G = 0;
		my @info = split(/\t+/, $pileup_line);
		my $chr = $info[0];
		my $loc = $info[1];
		$position = $chr."_$loc";
		if($info[3] >= 20){
			if($pileup_line =~ m/\*/ or $pileup_line =~ m/\+/ or $pileup_line =~ m/\-/){
				$mpileup_positions{"$position"} = "N";
			}
			else{
				for($genotype = 0; $genotype < length $info[4]; $genotype++){
					$allele = substr($info[4], $genotype, 1);
					if($allele eq "\." or $allele eq "\," or $allele eq "\>" or $allele eq "\<"){
						$reference++;
					}
					elsif($allele =~ "[Aa]"){
						$A++;
					}
					elsif($allele =~ "[Tt]"){
						$T++;
					}
					elsif($allele =~ "[Cc]"){
						$C++;
					}
					elsif($allele =~ "[Gg]"){
						$G++;
					}
				}
				$RAF = $reference/($reference+$A+$T+$C+$G);
				if($RAF >= 0.5){
					$mpileup_positions{"$position"} = $info[2];
				}
				elsif(($A > $T) && ($A > $C) && ($A > $G)){
					$mpileup_positions{"$position"} = "A";
				}
				elsif(($T > $A) && ($T > $C) && ($T > $G)){
					$mpileup_positions{"$position"} = "T";
				}
				elsif(($G > $A) && ($G > $C) && ($G > $T)){
					$mpileup_positions{"$position"} = "G";
				}
				elsif(($C > $A) && ($C > $G) && ($C > $T)){
					$mpileup_positions{"$position"} = "C";
				}
				else{$mpileup_positions{"$position"} = "N"};
			}
		}
		else{$mpileup_positions{"$position"} = "N"};
#		print "$mpileup_positions{$position}\n";
#		sleep(1);
	}
	foreach my $key (keys %var_pos){ 
    	if(defined $mpileup_positions{"$key"}){
    		print OUT "$mpileup_positions{$key}";
    	}
    	else{print OUT "N"}
	} 
	print OUT "\n";
}
close(IN);
close(OUT);

exit;