#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

mkdir "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_intergenic_files";

$ingen="${base_dir}/Data/Alignment_files/${species}_alignment.fasta";

$incoor="${base_dir}/Analysis/Gene_intergenic_coordinates/${species}_Gene_intergenic_coordinates/${species}_intergenic_coordinates.tab";

$count=0;
open INGEN, $ingen;
while(<INGEN>){
	$line=$_;
	chomp $line;
	if($line =~ /^>(\S+)/){
		$id=$1;
	}elsif($line =~ /^([ATGCN]+)/){
		$seq=$1;
		@seq_array=split(//, $seq);
		open INCOOR, $incoor;
		while(<INCOOR>){
			$line=$_;
			chomp $line;
			if($line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+\d+\s+\S+/){
				$intergenic_id=$1;$fir=$2;$sec=$3;
				if($line !~ /^Name\tStart\tEnd\tLength\tType/){
					if($sec-$fir>=0){
						$intergenic="";
						$indfir=($fir-1);
						$indsec=($sec-1);
						for $x($indfir..$indsec){
							$intergenic=$intergenic.$seq_array[$x];
						}
					
						open OUTPUT, ">>${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_intergenic_files/$intergenic_id.fasta";
						print OUTPUT ">$id\n$intergenic\n";
					}
				}
			}
		}
		$count++;
		#print "Isolate $count completed.\n";
	}
}

print "$species intergenic alignment split.\n";

print LOG "$species intergenic alignment split.\n";

