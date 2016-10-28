#!/usr/bin/perl -w

use List::Util qw(shuffle);

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];
$threshold=$ARGV[3];

open LOG, ">>${base_dir}/Analysis/log.txt";

if($threshold == 95){
	$threshold_folder="";
}else{
	$threshold_folder="/threshold_$threshold";
}

open OUTPUT_A, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}$threshold_folder/${species}_core_gene_alignment_no_STO_half_a.fasta";
open OUTPUT_B, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}$threshold_folder/${species}_core_gene_alignment_no_STO_half_b.fasta";

open INPUT, "${base_dir}/Analysis/${analysis}/${species}_${analysis}$threshold_folder/${species}_core_gene_alignment_no_STO.fasta";
while(<INPUT>){
	if(/^([ATGCN]+)/){
		$seq=$1;
		$len=length($seq);
		last;
	}
}

$no_codons=($len/3);
$mid_point=int($no_codons/2);

#print "$no_codons\n";

@codon_array=();
for($i=0; $i<$no_codons; $i++){
	$codon_array[$i]=$i;
}

@shuffled_codon_array=shuffle(@codon_array);

$count=0;
open INPUT, "${base_dir}/Analysis/${analysis}/${species}_${analysis}$threshold_folder/${species}_core_gene_alignment_no_STO.fasta";
while(<INPUT>){
	if(/^>(\S+)/){
		$id=$1;
		print OUTPUT_A ">$id\n";
		print OUTPUT_B ">$id\n";
	}elsif(/^([ATGCN]+)/){
		$seq=$1;
		@seq_array=split(//, $seq);
		
		for($i=0; $i<$mid_point; $i++){
			$pos=($shuffled_codon_array[$i]*3);
			$codon="$seq_array[$pos]$seq_array[($pos+1)]$seq_array[($pos+2)]";
			
			print OUTPUT_A "$codon";
		}
		print OUTPUT_A "\n";
		
		for($i=$mid_point; $i<$no_codons; $i++){
			$pos=($shuffled_codon_array[$i]*3);
			$codon="$seq_array[$pos]$seq_array[($pos+1)]$seq_array[($pos+2)]";
			
			print OUTPUT_B "$codon";
		}
		print OUTPUT_B "\n";
		
		$count++;
		#print "Isolate $count completed.\n";
	}
}

print "$species core genes shuffled.\n";

print LOG "$species core genes shuffled.\n";

