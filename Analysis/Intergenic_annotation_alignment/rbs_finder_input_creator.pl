#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

#open LOG, ">>${base_dir}/Analysis/log.txt";

$incoor="${base_dir}/Analysis/Gene_intergenic_coordinates/${species}_Gene_intergenic_coordinates/${species}_gene_coordinates.tab";

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_gene_coordinates_for_rbs_finder.tab";

open INCOOR, $incoor;
while(<INCOOR>){
	$line=$_;
	chomp $line;
	if($line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/){
		$gene_id=$1;$fir=$2;$sec=$3;$len=$4;$type=$5;$dir=$6;
		if($line !~ /^Name\tStart\tEnd\tLength\tType/){
			if($sec-$fir>=0){
				if($len%3==0 && $type eq "CDS"){
					if($dir eq "Reverse"){
						print OUTPUT "$gene_id\t$sec\t$fir\n";
					}else{
						print OUTPUT "$gene_id\t$fir\t$sec\n";
					}
				}
			}
		}
	}
}

