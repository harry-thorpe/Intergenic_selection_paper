#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_promoter_match.tab";

open IN_COOR, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_core_intergenics.tab";
while(<IN_COOR>){
	$line=$_;
	chomp $line;
	if($line =~ /^\S+\s+(\d+)\s+(\d+)\s+\d+\s+\S+/){
		$sta=$1;$end=$2;
		if($line !~ /^Name\tStart\tEnd\tLength\tType/){
			for($i=$sta; $i<=$end; $i++){
				$intergenic_hash{$i}=1;
				$intergenic_sta_hash{$i}=($sta - 1);
				$intergenic_end_hash{$i}=($end + 1);
			}
		}
	}
}

open IN_PRO, "${base_dir}/Data/Promoter_files/${species}_promoters.tab";
while(<IN_PRO>){
	$line=$_;
	chomp $line;
	@line_array=split(/\t/, $line, -1);
	
	if($line !~ /^ID\tSeqName/){
		
		$sta=$line_array[3];
		$end=$line_array[4];
		$strand=$line_array[7];
		# Start position is not correct in the file!
		$sta=($sta+1);
		
		for($i=$sta; $i<($sta + 6); $i++){
			$pro_hash{$i}=$strand;
		}
		
		for($i=$end; $i>($end - 6); $i--){
			$pro_hash{$i}=$strand;
		}
	}
}

open INCOOR, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_core_genes.tab";
while(<INCOOR>){
	$line=$_;
	chomp $line;
	if($line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/){
		$gene_id=$1;$fir=$2;$sec=$3;$len=$4;$type=$5;$dir=$6;
		if($line !~ /^Name\tStart\tEnd\tLength\tType/){
			if($sec-$fir>=0){
				if($len%3==0 && $type eq "CDS"){
					for $x($fir..$sec){
						$gene_hash{$x}=$dir;
						$gene_name_hash{$x}=$gene_id;
					}
				}
			}
		}
	}
}

open INPUT, "${base_dir}/Data/GFF_files/$species.gff";
while(<INPUT>){
	$line=$_;
	chomp $line;
	@line_array=split(/\t/, $line, -1);
	
	$sta="";
	$end="";
	if($line !~ /^##/){
		if($line_array[2] ne "sig_peptide" && $line_array[2] ne "misc_RNA"){
			$sta=$line_array[3];
			$end=$line_array[4];
			
			for $x($sta..$end){
				$gene_info_hash{$x}=$line;
			}
		}
	}elsif($line =~ /^##FASTA/){
		last;
	}
}

#isolates

open INPUT, "${base_dir}/Analysis/Mutation_intergenic_annotation/${species}_Mutation_intergenic_annotation/${species}_SNPs_intergenic.tab";
while(<INPUT>){
	$line=$_;
	chomp $line;
	@line_array=split(/\t/, $line, -1);
	
	$pos=$line_array[0];
	
	$gff_gene_info="";
	$gff_gene_product="";
	
	if($line !~ /^Position\tCategory/){
		if($line_array[1] eq "Promoter"){
			if($intergenic_hash{$pos}){
				if($pro_hash{$pos}){
					if($pro_hash{$pos} eq "+"){
						if($gene_hash{$intergenic_end_hash{$pos}}){
							if($gene_hash{$intergenic_end_hash{$pos}} eq "Forward"){
								if($gene_info_hash{$intergenic_end_hash{$pos}} =~ /gene=([^;]+)/){
									$gff_gene_info=$1;
								}
								if($gene_info_hash{$intergenic_end_hash{$pos}} =~ /product=(.+)/){
									$gff_gene_product=$1;
									
									$gff_gene_product=~tr/ /_/;
									$gff_gene_product=~tr/\./_/;
									$gff_gene_product=~tr/,/_/;
								}
								
								print OUTPUT "$pos\t$line_array[5]\t$gene_name_hash{$intergenic_end_hash{$pos}}\t$gff_gene_info\t$gff_gene_product\n";
							}else{
								#print OUTPUT "$pos\tGene_direction_conflict\n";
							}
						}else{
							#print OUTPUT "$pos\tGene_not_present\n";
						}
					}elsif($pro_hash{$pos} eq "-"){
						if($gene_hash{$intergenic_sta_hash{$pos}}){
							if($gene_hash{$intergenic_sta_hash{$pos}} eq "Reverse"){
								if($gene_info_hash{$intergenic_sta_hash{$pos}} =~ /gene=([^;]+)/){
									$gff_gene_info=$1;
								}
								if($gene_info_hash{$intergenic_sta_hash{$pos}} =~ /product=(.+)/){
									$gff_gene_product=$1;
									
									$gff_gene_product=~tr/ /_/;
									$gff_gene_product=~tr/\./_/;
									$gff_gene_product=~tr/,/_/;
								}
								
								print OUTPUT "$pos\t$line_array[5]\t$gene_name_hash{$intergenic_sta_hash{$pos}}\t$gff_gene_info\t$gff_gene_product\n";
							}else{
								#print OUTPUT "$pos\tGene_direction_conflict\n";
							}
						}else{
							#print OUTPUT "$pos\tGene_not_present\n";
						}
					}
				}else{
					print OUTPUT "$pos\tNot_promoter\n";
				}
			}else{
				print OUTPUT "$pos\tNot_intergenic\n";
			}
		}
	}
}

print "$species promoter genes matched.\n";

print LOG "$species promoter genes matched.\n";

