#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_dnds_simulated.csv";

print OUTPUT "Isolate_1,Isolate_2,dN,dS,dN/dS\n";

$include=0;
$dnds="";
open INPUT, "$base_dir/Analysis/$analysis/${species}_${analysis}/dnds_tmp/dnds_out.txt";
while(<INPUT>){
	$line=$_;
	chomp $line;
	
	if($line =~ /^Use runmode = -2 for ML pairwise comparison.\)/){
		$include=1;
	}elsif($line =~ /^\(B\) Yang & Nielsen \(2000\) method/){
		$include=0;
		last;
	}elsif($line =~ /\S+/){
		if($include == 1){
			if($line =~ /^(\S+)\s+(\S+.*)/){
				$id=$1;$values=$2;
				
				push @id_array, $id;
				
				$values=~s/\) / /g;
				$values=~s/\(//g;
				$values=~s/\)//g;
				
				if($dnds eq ""){
					$dnds=$values;
				}else{
					$dnds="$dnds $values";
				}
				
			}elsif($line =~ /^(\S+)/){
				$id=$1;
				
				push @id_array, $id;
			}
		}
	}
}

$id_count=scalar(@id_array);

@dnds_array=split(/\s+/, $dnds);

$dnds_count=scalar(@dnds_array);

if((($id_count * ($id_count - 1)) / 2) == ($dnds_count / 3)){
	print "all good\n";
	
	$count=0;
	for($i=0; $i<$id_count; $i++){
		for($j=0; $j<$id_count; $j++){
			if($i > $j){
				$dnds=$dnds_array[($count*3)];
				$dn=$dnds_array[(($count*3)+1)];
				$ds=$dnds_array[(($count*3)+2)];
				
				print OUTPUT "$id_array[$i],$id_array[$j],$dn,$ds,$dnds\n";
				
				$count++;
			}
		}
	}
}

