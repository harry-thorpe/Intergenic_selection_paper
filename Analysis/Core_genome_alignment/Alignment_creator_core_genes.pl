#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

$ingen="${base_dir}/Data/Alignment_files/${species}_alignment.fasta";

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_core_gene_alignment.fasta";

$count=0;
open INGEN, $ingen;
while(<INGEN>){
	$line=$_;
	chomp $line;
	if($line =~ /^>(\S+)/){
		$id=$1;
		
		print OUTPUT ">$id\n";
	}elsif($line =~ /^([ATGCN]+)/){
		$seq=$1;
		@seq_array=split(//, $seq);
		open INCOOR, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_core_genes.tab";
		while(<INCOOR>){
			$line=$_;
			chomp $line;
			if($line =~ /^\S+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/){
				$fir=$1;$sec=$2;$len=$3;$type=$4;$dir=$5;
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
						
							print OUTPUT "$gene";
						}
					}
				}
			}
		}
		
		print OUTPUT "\n";
		
		$count++;
		#print "Isolate $count completed.\n";
	}
}

print "$species core gene alignment created.\n";

print LOG "$species core gene alignment created.\n";

