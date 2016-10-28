#!/usr/bin/perl -w

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

mkdir "${base_dir}/Analysis/${analysis}/${species}_${analysis}$threshold_folder/Gene_files_half_a";
mkdir "${base_dir}/Analysis/${analysis}/${species}_${analysis}$threshold_folder/Gene_files_half_b";

$count=0;
open INPUT, "${base_dir}/Analysis/${analysis}/${species}_${analysis}$threshold_folder/${species}_core_gene_alignment_no_STO_half_a.fasta";
while(<INPUT>){
	$line=$_;
	chomp $line;
	
	if($line =~/^>(\S+)/){
		$id=$1;
	}elsif($line =~ /^([ATGCN]+)/){
		$seq=$1;
		
		open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}$threshold_folder/Gene_files_half_a/$id.fasta";
		print OUTPUT ">$id\n$seq\n";
		
		$count++;
		#print "Isolate $count completed.\n";
	}
}

$count=0;
open INPUT, "${base_dir}/Analysis/${analysis}/${species}_${analysis}$threshold_folder/${species}_core_gene_alignment_no_STO_half_b.fasta";
while(<INPUT>){
	$line=$_;
	chomp $line;
	
	if($line =~/^>(\S+)/){
		$id=$1;
	}elsif($line =~ /^([ATGCN]+)/){
		$seq=$1;
		
		open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}$threshold_folder/Gene_files_half_b/$id.fasta";
		print OUTPUT ">$id\n$seq\n";
		
		$count++;
		#print "Isolate $count completed.\n";
	}
}

print "$species core gene alignment split.\n";

print LOG "$species core gene alignment split.\n";

