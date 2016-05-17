#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

mkdir "${base_dir}/Analysis/${analysis}/${species}_${analysis}/Intergenic_files";

$count=0;
open INPUT, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_core_intergenic_alignment.fasta";
while(<INPUT>){
	$line=$_;
	chomp $line;
	
	if($line =~/^>(\S+)/){
		$id=$1;
	}elsif($line =~ /^([ATGCN]+)/){
		$seq=$1;
		
		open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/Intergenic_files/$id.fasta";
		print OUTPUT ">$id\n$seq\n";
		
		$count++;
		#print "Isolate $count completed.\n";
	}
}

print "$species core intergenic alignment split.\n";

print LOG "$species core intergenic alignment split.\n";

