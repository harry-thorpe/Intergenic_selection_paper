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

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}$threshold_folder/${species}_PSM.csv";

print OUTPUT "Category,Total_SNPs,Singletons,Doubletons,Proportion_of_Singletons,Proportion_of_Doubletons\n";

open INPUT, "${base_dir}/Analysis/${analysis}/${species}_${analysis}$threshold_folder/${species}_mutations.tab";
while(<INPUT>){
	$line=$_;
	chomp $line;
	@line_array=split(/\t/, $line);
	
	if($line !~ /^Position\tCategory/){
		$category=$line_array[1];
		$freq=$line_array[5];
		
		$mut_type_hash{$category}++;
		
		if($freq == 1){
			$mut_type_singleton_hash{$category}++;
		}elsif($freq == 2){
			$mut_type_doubleton_hash{$category}++;
		}
	}
}

@category_array=("Synonymous", "Nonsynonymous", "Nonsense", "Intergenic");

foreach $category(@category_array){
	$PSM=($mut_type_singleton_hash{$category} / $mut_type_hash{$category});
	$PDM=($mut_type_doubleton_hash{$category} / ($mut_type_hash{$category} - $mut_type_singleton_hash{$category}));
	
	print OUTPUT "$category,$mut_type_hash{$category},$mut_type_singleton_hash{$category},$mut_type_doubleton_hash{$category},$PSM,$PDM\n";
}

print "$species PSM analysis completed.\n";

print LOG "$species PSM analysis completed.\n";

