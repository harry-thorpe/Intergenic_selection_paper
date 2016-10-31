#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_PSM_intergenic.csv";

print OUTPUT "Category,Total_SNPs,Singletons,Proportion_of_Singletons\n";

open OUTPUT_2, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_SNPs_intergenic.tab";

print OUTPUT_2 "Position\tCategory\tReference\tSNP\tReference_count\tSNP_count\tReference_isolates\tSNP_isolates\n";


open IN_COOR, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_core_intergenics.tab";
while(<IN_COOR>){
	$line=$_;
	chomp $line;
	if($line =~ /^\S+\s+(\d+)\s+(\d+)\s+\d+\s+\S+/){
		$sta=$1;$end=$2;
		$sta=($sta-1);
		$end=($end-1);
		if($line !~ /^Name\tStart\tEnd\tLength\tType/){
			for($i=$sta; $i<=$end; $i++){
				$intergenic_hash{$i}=1;
				
				push @intergenic_array, $i;
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
		$sta=($sta-1);
		$end=($end-1);
		# Start position is not correct in the file!
		$sta=($sta+1);
		
		for($i=$sta; $i<($sta + 6); $i++){
			$pro_hash{$i}=1;
		}
		
		for($i=$end; $i>($end - 6); $i--){
			$pro_hash{$i}=1;
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
		$sta=($sta-1);
		$end=($end-1);
		
		for($i=$sta; $i<=$end; $i++){
			$ter_hash{$i}=1;
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
			$sta=($sta-1);
			$end=($end-1);
		
			for($i=$sta; $i<=$end; $i++){
				$rna_hash{$i}=1;
			}
		}
	}elsif($line =~ /^##FASTA/){
		last;
	}
}

$count=0;
open INGEN, "${base_dir}/Data/Alignment_files/${species}_alignment.fasta";
while(<INGEN>){
	$line=$_;
	chomp $line;
	
	if($line =~ /^>(\S+)/){
		$id=$1;
		
	}elsif($line =~ /^([ATGCN]+)/){
	
		@seq_array=split(//, $line);
	
		foreach $i(@intergenic_array){
			if($seq_array[$i] !~ /N/){
				if($intergenic_hash{$i}){
					
					$base_hash{$i}{$seq_array[$i]}++;
				}
			}
		}
	
		$count++;
		#print "Isolate $count completed.\n";
	}
}

@base_pos_array=keys(%base_hash);
@base_pos_array=sort { $a <=> $b } @base_pos_array;

@biallelic_sites_array=();
foreach $base_pos(@base_pos_array){
	@base_array=sort { $base_hash{$base_pos}{$b} <=> $base_hash{$base_pos}{$a} } keys(%{$base_hash{$base_pos}});

	$base_count=scalar(@base_array);

	if($base_count == 2){
		push @biallelic_sites_array, $base_pos;
	}
}

%base_hash=();

$count=0;
open INGEN, "${base_dir}/Data/Alignment_files/${species}_alignment.fasta";
while(<INGEN>){
	$line=$_;
	chomp $line;
	
	if($line =~ /^>(\S+)/){
		$id=$1;
		
	}elsif($line =~ /^([ATGCN]+)/){
	
		@seq_array=split(//, $line);
	
		foreach $i(@biallelic_sites_array){
			if($seq_array[$i] !~ /N/){
				if($intergenic_hash{$i} && !$pro_hash{$i} && !$ter_hash{$i} && !$rna_hash{$i}){
					
					$mut="Unannotated";
					
					$base_hash{$mut}{$i}{$seq_array[$i]}++;
					
					if(!$base_isolate_hash{$mut}{$i}{$seq_array[$i]}){
						$base_isolate_hash{$mut}{$i}{$seq_array[$i]}=$id;
					}else{
						$base_isolate_hash{$mut}{$i}{$seq_array[$i]}="$base_isolate_hash{$mut}{$i}{$seq_array[$i]},$id";
					}
				}
				
				if($intergenic_hash{$i} && $pro_hash{$i}){
					
					$mut="Promoter";
					
					$base_hash{$mut}{$i}{$seq_array[$i]}++;
					
					if(!$base_isolate_hash{$mut}{$i}{$seq_array[$i]}){
						$base_isolate_hash{$mut}{$i}{$seq_array[$i]}=$id;
					}else{
						$base_isolate_hash{$mut}{$i}{$seq_array[$i]}="$base_isolate_hash{$mut}{$i}{$seq_array[$i]},$id";
					}
				}
				
				if($intergenic_hash{$i} && $ter_hash{$i}){
					
					$mut="Terminator";
					
					$base_hash{$mut}{$i}{$seq_array[$i]}++;
					
					if(!$base_isolate_hash{$mut}{$i}{$seq_array[$i]}){
						$base_isolate_hash{$mut}{$i}{$seq_array[$i]}=$id;
					}else{
						$base_isolate_hash{$mut}{$i}{$seq_array[$i]}="$base_isolate_hash{$mut}{$i}{$seq_array[$i]},$id";
					}
				}
				
				if($intergenic_hash{$i} && $rna_hash{$i}){
					
					$mut="Non_coding_RNA";
					
					$base_hash{$mut}{$i}{$seq_array[$i]}++;
					
					if(!$base_isolate_hash{$mut}{$i}{$seq_array[$i]}){
						$base_isolate_hash{$mut}{$i}{$seq_array[$i]}=$id;
					}else{
						$base_isolate_hash{$mut}{$i}{$seq_array[$i]}="$base_isolate_hash{$mut}{$i}{$seq_array[$i]},$id";
					}
				}
			}
		}
	
		$count++;
		#print "Isolate $count completed.\n";
	}
}

@category_array=("Non_coding_RNA", "Promoter", "Terminator", "Unannotated");

foreach $category(@category_array){

	@base_pos_array=keys(%{$base_hash{$category}});
	@base_pos_array=sort { $a <=> $b } @base_pos_array;

	foreach $base_pos(@base_pos_array){
		@base_array=sort { $base_hash{$category}{$base_pos}{$b} <=> $base_hash{$category}{$base_pos}{$a} } keys(%{$base_hash{$category}{$base_pos}});

		$base_count=scalar(@base_array);

		if($base_count == 2){
	
			$mut="$category";
	
			$mut_type_hash{$mut}++;
	
			if($base_hash{$category}{$base_pos}{$base_array[1]} == 1){
				$mut_type_singleton_hash{$mut}++;
			}
			
			$real_base_pos=($base_pos + 1);
		
			print OUTPUT_2 "$real_base_pos\t$mut\t$base_array[0]\t$base_array[1]\t$base_hash{$category}{$base_pos}{$base_array[0]}\t$base_hash{$category}{$base_pos}{$base_array[1]}\t$base_isolate_hash{$category}{$base_pos}{$base_array[0]}\t$base_isolate_hash{$category}{$base_pos}{$base_array[1]}\n";
		
		}
	}
}

@mut_array=keys(%mut_type_hash);
@mut_array=sort(@mut_array);

foreach $mut(@mut_array){

	$prop=($mut_type_singleton_hash{$mut} / $mut_type_hash{$mut});

	print OUTPUT "$mut,$mut_type_hash{$mut},$mut_type_singleton_hash{$mut},$prop\n";
}

print "$species PSM intergenic annotation analysis completed.\n";

print LOG "$species PSM intergenic annotation analysis completed.\n";

