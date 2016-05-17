#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

mkdir "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_gene_files";

$ingen="${base_dir}/Data/Alignment_files/${species}_alignment.fasta";

$incoor="${base_dir}/Analysis/Gene_intergenic_coordinates/${species}_Gene_intergenic_coordinates/${species}_gene_coordinates.tab";

open OUTPUT_ISOLATE, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_isolates.txt";

$count=0;
open INGEN, $ingen;
while(<INGEN>){
	$line=$_;
	chomp $line;
	if($line =~ /^>(\S+)/){
		$id=$1;
		
		print OUTPUT_ISOLATE "$id\n";
	}elsif($line =~ /^([ATGCN]+)/){
		$seq=$1;
		@seq_array=split(//, $seq);
		open INCOOR, $incoor;
		while(<INCOOR>){
			$line=$_;
			chomp $line;
			if($line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/){
				$gene_id=$1;$fir=$2;$sec=$3;$len=$4;$type=$5;$dir=$6;
				if($line !~ /^Name\tStart\tEnd\tLength\tType/){
					if($sec-$fir>=0){
						if($len%3==0 && $type eq "CDS"){
							$gene="";
							$indfir=($fir-1);
							$indsec=($sec-1);
							for $x($indfir..$indsec){
								$gene=$gene.$seq_array[$x];
							}
							if($dir eq "Reverse"){
								$gene=~tr/ATGC/TACG/;
								$gene=reverse($gene);
							}
						
							open OUTPUT, ">>${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_gene_files/$gene_id.fasta";
							print OUTPUT ">$id\n$gene\n";
						}
					}
				}
			}
		}
		$count++;
		#print "Isolate $count completed.\n";
	}
}

print "$species gene alignment split.\n";

print LOG "$species gene alignment split.\n";

