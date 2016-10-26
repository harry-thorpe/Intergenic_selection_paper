#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_variants_summarised.tab";

print OUTPUT "Variant\tScore\tDepth\tAF1\tMQ\tSupport\n";

$count=0;
open INPUT, "${base_dir}/Data/Alignment_files/${species}_alignment.fasta";
while(<INPUT>){
	$line=$_;
	chomp $line;
	
	if($line =~ /^>(\S+)/){
		$id=$1;
		
		if($count > 0){
			push @id_array, $id;
		}
		
		$count++;
	}
}
$id_count=scalar(@id_array);

for($i=0; $i<$id_count; $i++){
	$id=$id_array[$i];
	
	open INPUT, "${base_dir}/Data/vcf_files/${species}_vcf_files/${id}.variants.vcf";
	while(<INPUT>){
		$line=$_;
		chomp $line;
		
		if($line !~ /^#/ && $line !~ /INDEL/){
			if($line =~ /^\S+\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)/){
				$pos=$1;$ref=$2;$alt=$3;$score=$4;
			}
			if($line =~ /DP=(\d+)/){
				$dp=$1;
			}
			if($line =~ /AF1=(\d+)/){
				$af1=$1;
			}
			if($line =~ /MQ=(\d+)/){
				$mq=$1;
			}
			if($line =~ /DP4=\d+,\d+,(\d+),(\d+)/){
				$alt1=$1;$alt2=$2;
			}
			$support=(($alt1+$alt2)/$dp);
		
			print OUTPUT "${pos}_${ref}_${alt}\t$score\t$dp\t$af1\t$mq\t$support\n";
		}
	}
	#print "$i $id_array[$i] completed.\n";
}
close OUTPUT;

open INPUT, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_variants_summarised.tab";
while(<INPUT>){
	$line=$_;
	chomp $line;
	@line_array=split(/\t/, $line);
	
	if($line !~ /^Variant/){
		if(!$variant_hash{$line_array[0]}){
			$variant_hash{$line_array[0]}=1;
		}else{
			$variant_hash{$line_array[0]}++;
		}
	}
}

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_variants_summarised_counted.tab";

print OUTPUT "Variant\tScore\tDepth\tAF1\tMQ\tSupport\tCount\n";

open INPUT, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_variants_summarised.tab";
while(<INPUT>){
	$line=$_;
	chomp $line;
	@line_array=split(/\t/, $line);
	
	if($line !~ /^Variant/){
		print OUTPUT "$line\t$variant_hash{$line_array[0]}\n";
	}
}
close OUTPUT;

