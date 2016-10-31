#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

open OUTPUT_RBS, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_core_intergenic_rbs_alignment.fasta";
open OUTPUT_TER, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_core_intergenic_terminator_alignment.fasta";
open OUTPUT_PRO, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_core_intergenic_promoter_alignment.fasta";
open OUTPUT_RNA, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_core_intergenic_non_coding_RNA_alignment.fasta";
open OUTPUT_UNA, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_core_intergenic_unannotated_alignment.fasta";

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
			}
		}
	}
}

open IN_RBS, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_rbs_coordinates.tab";
while(<IN_RBS>){
	$line=$_;
	chomp $line;
	
	if($line =~ /^\S+\s+(\d+)\s+(\d+)/){
		$sta=$1;
		$end=$2;
		$sta=($sta-1);
		$end=($end-1);
		
		for($i=$sta; $i<=$end; $i++){
			$rbs_hash{$i}=1;
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
		
		print OUTPUT_RBS ">$id\n";
		print OUTPUT_TER ">$id\n";
		print OUTPUT_PRO ">$id\n";
		print OUTPUT_RNA ">$id\n";
		print OUTPUT_UNA ">$id\n";
	}elsif($line =~ /^([ATGCN]+)/){
		$seq=$1;
		@seq_array=split(//, $seq);
		$genome_len=length($seq);
		
		for($i=0; $i<$genome_len; $i++){
			if($intergenic_hash{$i} && $rbs_hash{$i}){
				print OUTPUT_RBS "$seq_array[$i]";
			}elsif($intergenic_hash{$i} && $ter_hash{$i}){
				print OUTPUT_TER "$seq_array[$i]";
			}elsif($intergenic_hash{$i} && $pro_hash{$i}){
				print OUTPUT_PRO "$seq_array[$i]";
			}elsif($intergenic_hash{$i} && $rna_hash{$i}){
				print OUTPUT_RNA "$seq_array[$i]";
			}elsif($intergenic_hash{$i} && !$rbs_hash{$i} && !$ter_hash{$i} && !$pro_hash{$i} && !$rna_hash{$i}){
				print OUTPUT_UNA "$seq_array[$i]";
			}
		}
		
		print OUTPUT_RBS "\n";
		print OUTPUT_TER "\n";
		print OUTPUT_PRO "\n";
		print OUTPUT_RNA "\n";
		print OUTPUT_UNA "\n";
		
		$count++;
		#print "Isolate $count completed.\n";
	}
}

print "$species intergenic annotation alignments created.\n";

print LOG "$species intergenic annotation alignments created.\n";

