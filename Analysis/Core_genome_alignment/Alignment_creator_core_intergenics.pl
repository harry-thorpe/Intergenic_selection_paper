#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];
$threshold=$ARGV[3];

open LOG, ">>${base_dir}/Analysis/log.txt";

$ingen="${base_dir}/Data/Alignment_files/${species}_alignment.fasta";

if($threshold == 95){
	$threshold_folder="";
}else{
	$threshold_folder="/threshold_$threshold";
}

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}$threshold_folder/${species}_core_intergenic_alignment.fasta";

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
		open INCOOR, "${base_dir}/Analysis/${analysis}/${species}_${analysis}$threshold_folder/${species}_core_intergenics.tab";
		while(<INCOOR>){
			$line=$_;
			chomp $line;
			if($line =~ /^\S+\s+(\d+)\s+(\d+)\s+\d+\s+\S+/){
				$fir=$1;$sec=$2;
				if($line !~ /^Name\tStart\tEnd\tLength\tType/){
					if($sec-$fir>=0){
						$intergenic="";
						$indfir=($fir-1);
						$indsec=($sec-1);
						for $x($indfir..$indsec){
							$intergenic=$intergenic.$seq_array[$x];
						}
				
						print OUTPUT "$intergenic";
					}
				}
			}
		}
	
		print OUTPUT "\n";
	
		$count++;
		#print "Isolate $count completed.\n";
	}
}

print "$species core intergenic alignment created.\n";

print LOG "$species core intergenic alignment created.\n";

