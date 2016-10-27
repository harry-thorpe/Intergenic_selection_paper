#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_mutation_intergenic_unannotated_distance.tab";

print OUTPUT "Category\tDistance\tSNPs\tSites\tSNP_density\n";

$step_size=10;
$max=150;

open IN_PRO, "${base_dir}/Data/Promoter_files/${species}_promoters.tab";
while(<IN_PRO>){
	$line=$_;
	chomp $line;
	@line_array=split(/\t/, $line, -1);
	
	if($line !~ /^ID\tSeqName/){
		
		$sta=$line_array[3];
		$end=$line_array[4];
		# Start position is not correct in the file!
		$sta=($sta+1);
		
		for($i=$sta; $i<($sta + 6); $i++){
			$mask_sites_hash{$i}=1;
		}
		
		for($i=$end; $i>($end - 6); $i--){
			$mask_sites_hash{$i}=1;
		}
	}
}

open IN_TER, "${base_dir}/Data/Terminator_files/${species}_terminators.tab";
while(<IN_TER>){
	$line=$_;
	chomp $line;
	
	if($line =~ /\S+\s+(\d+)\s+\.\.\s+(\d+)/){
		$sta=$1;
		$end=$2;
		
		for($i=$sta; $i<=$end; $i++){
			$mask_sites_hash{$i}=1;
		}
	}
}

open IN_RNA, "${base_dir}/Data/GFF_files/$species.gff";
while(<IN_RNA>){
	$line=$_;
	chomp $line;
	@line_array=split(/\t/, $line, -1);
	
	$sta="";
	$end="";
	
	if($line !~ /^##/){
		if($line_array[2] eq "misc_RNA"){
			$sta=$line_array[3];
			$end=$line_array[4];
		
			for($i=$sta; $i<=$end; $i++){
				$mask_sites_hash{$i}=1;
			}
		}
	}elsif($line =~ /^##FASTA/){
		last;
	}
}

open IN_COOR, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_core_intergenics.tab";
while(<IN_COOR>){
	$line=$_;
	chomp $line;
	if($line =~ /^\S+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)/){
		$sta=$1;$end=$2;$real_len=$3;$type=$4;
		if($line !~ /^Name\tStart\tEnd\tLength\tType/){
			for($i=$sta; $i<=$end; $i++){
				$int_sta_hash{$i}=$sta;
				$int_end_hash{$i}=$end;
				$int_typ_hash{$i}=$type;
			}
			
			if($type eq "Intergenic_co-oriented_F"){
				$len=((int($real_len/10))*10);
		
				for($i=0; $i<=$len; $i+=$step_size){
					for($j=$i; $j<($i+$step_size); $j++){
						if(!$mask_sites_hash{$j}){
							$three_dis_site_hash{$i}++;
						}
					}
			
					for($j=$i; $j<($i+$step_size); $j++){
						if(!$mask_sites_hash{$j}){
							$five_dis_site_hash{($len-$i)}++;
						}
					}
				}
			}elsif($type eq "Intergenic_co-oriented_R"){
				$len=((int($real_len/10))*10);
		
				for($i=0; $i<=$len; $i+=$step_size){
					for($j=$i; $j<($i+$step_size); $j++){
						if(!$mask_sites_hash{$j}){
							$five_dis_site_hash{$i}++;
						}
					}
			
					for($j=$i; $j<($i+$step_size); $j++){
						if(!$mask_sites_hash{$j}){
							$three_dis_site_hash{($len-$i)}++;
						}
					}
				}
			}elsif($type eq "Intergenic_double_promoter"){
				$len=((int($real_len/10))*10);
		
				for($i=0; $i<=$len; $i+=10){
					#$five_dis_site_hash{$i}++;
					#$five_dis_site_hash{$i}++;
				}
			}elsif($type eq "Intergenic_double_terminator"){
				$len=((int($real_len/10))*10);
		
				for($i=0; $i<=$len; $i+=10){
					#$three_dis_site_hash{$i}++;
					#$three_dis_site_hash{$i}++;
				}
			}
		}
	}
}

open IN_COOR, "${base_dir}/Analysis/Mutation_intergenic_annotation/${species}_Mutation_intergenic_annotation/${species}_SNPs_intergenic.tab";
while(<IN_COOR>){
	$line=$_;
	chomp $line;
	@line_array=split(/\t/, $line, -1);
	
	if($line !~ /^Position\tCategory\tReference\tSNP/){
		$snp=$line_array[0];
	
		if(!$mask_sites_hash{$snp}){
			if($int_typ_hash{$snp} eq "Intergenic_co-oriented_F"){
				$five_dis=$int_end_hash{$snp}-$snp;
				$five_dis=((int($five_dis/10))*10);
			
				if($five_dis <= $max){
					$five_dis_snp_hash{$five_dis}++;
				}
			
				$three_dis=$snp-$int_sta_hash{$snp};
				$three_dis=((int($three_dis/10))*10);
			
				if($three_dis <= $max){
					$three_dis_snp_hash{$three_dis}++;
				}
			}elsif($int_typ_hash{$snp} eq "Intergenic_co-oriented_R"){
				$five_dis=$snp-$int_sta_hash{$snp};
				$five_dis=((int($five_dis/10))*10);
			
				if($five_dis <= $max){
					$five_dis_snp_hash{$five_dis}++;
				}
			
				$three_dis=$int_end_hash{$snp}-$snp;
				$three_dis=((int($three_dis/10))*10);
			
				if($three_dis <= $max){
					$three_dis_snp_hash{$three_dis}++;
				}
			}elsif($int_typ_hash{$snp} eq "Intergenic_double_promoter"){
				$five_dis=$snp-$int_sta_hash{$snp};
				$five_dis=((int($five_dis/10))*10);
			
				if($five_dis <= $max){
					#$five_dis_snp_hash{$five_dis}++;
				}
			
				$five_dis=$int_end_hash{$snp}-$snp;
				$five_dis=((int($five_dis/10))*10);
			
				if($five_dis <= $max){
					#$five_dis_snp_hash{$five_dis}++;
				}
			}elsif($int_typ_hash{$snp} eq "Intergenic_double_terminator"){
				$three_dis=$snp-$int_sta_hash{$snp};
				$three_dis=((int($three_dis/10))*10);
			
				if($three_dis <= $max){
					#$three_dis_snp_hash{$three_dis}++;
				}
			
				$three_dis=$int_end_hash{$snp}-$snp;
				$three_dis=((int($three_dis/10))*10);
			
				if($three_dis <= $max){
					#$three_dis_snp_hash{$three_dis}++;
				}
			}
		}
	}
}

for($i=0; $i<=$max; $i+=$step_size){
	if(!$five_dis_snp_hash{$i}){
		$five_dis_snp_hash{$i}=0;
	}
	
	if(!$three_dis_snp_hash{$i}){
		$three_dis_snp_hash{$i}=0;
	}
	
	$five_density=($five_dis_snp_hash{$i}/$five_dis_site_hash{$i});
	$three_density=($three_dis_snp_hash{$i}/$three_dis_site_hash{$i});
	
	print OUTPUT "Gene_start_5'\t$i\t$five_dis_snp_hash{$i}\t$five_dis_site_hash{$i}\t$five_density\n";
	print OUTPUT "Gene_end_3'\t$i\t$three_dis_snp_hash{$i}\t$three_dis_site_hash{$i}\t$three_density\n";
}

print "$species Mutation intergenic unannotated distance completed.\n";

print LOG "$species Mutation intergenic unannotated distance completed.\n";

