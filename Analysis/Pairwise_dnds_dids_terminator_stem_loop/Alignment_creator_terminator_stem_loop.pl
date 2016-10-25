#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open OUTPUT_TER_STEM, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_core_intergenic_terminator_stem_alignment.fasta";
open OUTPUT_TER_LOOP, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_core_intergenic_terminator_loop_alignment.fasta";

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

open IN_TER, "${base_dir}/Data/Terminator_files/${species}_terminators.tab";
while(<IN_TER>){
	$line=$_;
	chomp $line;
	
	if($line =~ /\S+\s+(\d+)\s+\.\.\s+(\d+)\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)/){
		$sta=$1;
		$end=$2;
		$sta=($sta-1);
		$end=($end-1);
		$lhs=$3;
		$loop=$4;
		$rhs=$5;
		
		$lhs=~s/-//g;
		$rhs=~s/-//g;
		$loop=~s/-//g;
		
		$lhs_len=length($lhs);
		$rhs_len=length($rhs);
		$loop_len=length($loop);
		
		$sta=$sta;
		$end=($sta+$lhs_len);
		for($i=$sta; $i<$end; $i++){
			$ter_stem_hash{$i}=1;
		}
		
		$sta=$end;
		$end=($end+$loop_len);
		for($i=$sta; $i<$end; $i++){
			$ter_loop_hash{$i}=1;
		}
		
		$sta=$end;
		$end=($sta+$rhs_len);
		for($i=$sta; $i<$end; $i++){
			$ter_stem_hash{$i}=1;
		}
	}
}

$count=0;
open INGEN, "${base_dir}/Data/Alignment_files/${species}_alignment.fasta";
while(<INGEN>){
	$line=$_;
	chomp $line;
	
	if($line =~ /^>(\S+)/){
		$id=$1;
		
		print OUTPUT_TER_STEM ">$id\n";
		print OUTPUT_TER_LOOP ">$id\n";
	}elsif($line =~ /^([ATGCN]+)/){
		$seq=$1;
		@seq_array=split(//, $seq);
		$genome_len=length($seq);
		
		for($i=0; $i<$genome_len; $i++){
			if($intergenic_hash{$i} && $ter_stem_hash{$i}){
				print OUTPUT_TER_STEM "$seq_array[$i]";
			}elsif($intergenic_hash{$i} && $ter_loop_hash{$i}){
				print OUTPUT_TER_LOOP "$seq_array[$i]";
			}
		}
		
		print OUTPUT_TER_STEM "\n";
		print OUTPUT_TER_LOOP "\n";
		
		$count++;
		#print "Isolate $count completed.\n";
	}
}

print "$species terminator stem loop alignments created.\n";

