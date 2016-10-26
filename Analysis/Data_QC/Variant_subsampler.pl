#!/usr/bin/perl -w

use List::Util 'shuffle';

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_variants_subsampled.tab";

$count=0;
@variant_array=();
open INPUT, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_variants_summarised_counted.tab";
while(<INPUT>){
	$line=$_;
	chomp $line;
	
	if($line !~ /^Variant/){
		$count++;
		
		push @variant_array, $count;
	}
}

$max=1000000;

if($count > $max){
	$step=int($count / $max);
}else{
	$step=1;
}

open INPUT, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_variants_summarised_counted.tab";
while(<INPUT>){
	$line=$_;
	chomp $line;
	
	if($line !~ /^Variant/){
		$count++;
		
		if($count % $step == 0){
			print OUTPUT "$line\n";
		}
	}else{
		print OUTPUT "$line\n";
	}
}

