#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];
$upstream_len=$ARGV[3];

open LOG, ">>${base_dir}/Analysis/log.txt";

mkdir "${base_dir}/Analysis/${analysis}/${species}_${analysis}/Intergenic_files_upstream_$upstream_len";

$count=0;
open INPUT, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_core_intergenic_alignment_upstream_$upstream_len.fasta";
while(<INPUT>){
	$line=$_;
	chomp $line;
	
	if($line =~/^>(\S+)/){
		$id=$1;
	}elsif($line =~ /^([ATGCN]+)/){
		$seq=$1;
		
		open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/Intergenic_files_upstream_$upstream_len/$id.fasta";
		print OUTPUT ">$id\n$seq\n";
		
		$count++;
		#print "Isolate $count completed.\n";
	}
}

print "$species core intergenic alignment upstream $upstream_len split.\n";

print LOG "$species core intergenic alignment upstream $upstream_len split.\n";

