#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];
$upstream_len=$ARGV[3];

open LOG, ">>${base_dir}/Analysis/log.txt";

$ingen="${base_dir}/Data/Alignment_files/${species}_alignment.fasta";

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_core_intergenic_alignment_upstream_$upstream_len.fasta";

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
		open INCOOR, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_core_intergenics.tab";
		while(<INCOOR>){
			$line=$_;
			chomp $line;
			if($line =~ /^\S+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)/){
				$sta=$1;$end=$2;$len=$3;$type=$4;
				if($line !~ /^Name\tStart\tEnd\tLength\tType/){
					if($end-$sta>=0){
						$intergenic="";
						$indsta=($sta-1);
						$indend=($end-1);
					
						if($type eq "Intergenic_co-oriented_F"){
							if($len >= $upstream_len){
								$new_sta=($indend - ($upstream_len - 1));
							
								for $x($new_sta..$indend){
									$intergenic=$intergenic.$seq_array[$x];
								}
							}else{
								for $x($indsta..$indend){
									$intergenic=$intergenic.$seq_array[$x];
								}
							}
						}elsif($type eq "Intergenic_co-oriented_R"){
							if($len >= $upstream_len){
								$new_end=($indsta + ($upstream_len - 1));
							
								for $x($indsta..$new_end){
									$intergenic=$intergenic.$seq_array[$x];
								}
							}else{
								for $x($indsta..$indend){
									$intergenic=$intergenic.$seq_array[$x];
								}
							}
						}elsif($type eq "Intergenic_double_promoter"){
							if($len >= ($upstream_len * 2)){
								$new_sta=($indend - ($upstream_len - 1));
							
								for $x($new_sta..$indend){
									$intergenic=$intergenic.$seq_array[$x];
								}
							
								$new_end=($indsta + ($upstream_len - 1));
							
								for $x($indsta..$new_end){
									$intergenic=$intergenic.$seq_array[$x];
								}
							}else{
								for $x($indsta..$indend){
									$intergenic=$intergenic.$seq_array[$x];
								}
							}
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

print "$species core intergenic alignment upstream $upstream_len created.\n";

print LOG "$species core intergenic alignment upstream $upstream_len created.\n";

