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

$threshold_d=($threshold/100);

%core_hash=();
@intergenic_array=();

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}$threshold_folder/${species}_core_intergenics.tab";
print OUTPUT "Name\tStart\tEnd\tLength\tType\n";

$incoor="${base_dir}/Analysis/Gene_intergenic_coordinates/${species}_Gene_intergenic_coordinates/${species}_intergenic_coordinates.tab";

open INCOOR, $incoor;
while(<INCOOR>){
	$line=$_;
	chomp $line;
	if($line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)/){
		$intergenic_id=$1;$fir=$2;$sec=$3;$len=$4;$type=$5;
		if($line !~ /^Name\tStart\tEnd\tLength\tType/){
			if($sec-$fir>=0){
				push @intergenic_array, $intergenic_id;
			}
		}
	}
}

foreach $intergenic(@intergenic_array){

	$isolate_prop=0;
	$isolate_count=0;
	$isolate_include=0;
	open INGEN, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_intergenic_files/$intergenic.fasta";
	while(<INGEN>){
		$line=$_;
		chomp $line;
		if($line =~ /^>\S+/){
		
			$isolate_count++;
		}elsif($line =~ /^([ATGCN]+)/){
			$seq=$1;
			@seq_array=split(//, $seq);
		
			$len=length($seq);
		
			$base_count=0;
			$base_prop=0;
			for($i=0; $i<$len; $i++){
				if($seq_array[$i] ne "N"){
					$base_count++;
				}
			}
			$base_prop=($base_count/$len);
		
			if($base_prop >= 0.9){
				$isolate_include++;
			}
		}
	}

	$isolate_prop=($isolate_include/$isolate_count);

	if($isolate_prop >= $threshold_d){
		$core_hash{$intergenic}=1;
	}

	#print "Intergenic region $intergenic completed.\n";
}

open INCOOR, $incoor;
while(<INCOOR>){
	$line=$_;
	chomp $line;
	if($line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)/){
		$intergenic_id=$1;$fir=$2;$sec=$3;$len=$4;$type=$5;
		if($line !~ /^Name\tStart\tEnd\tLength\tType/){
			if($sec-$fir>=0){
				if($core_hash{$intergenic_id}){
					print OUTPUT "$line\n";
				}
			}
		}
	}
}

print "$species core intergenics found.\n";

print LOG "$species core intergenics found.\n";

