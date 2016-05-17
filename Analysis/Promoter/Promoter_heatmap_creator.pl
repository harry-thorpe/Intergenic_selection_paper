#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_promoter_heatmap.tab";

open INPUT, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_promoter_match.tab";
while(<INPUT>){
	$line=$_;
	chomp $line;
	@line_array=split(/\t/, $line, -1);
	
	push @mutation_array, $line_array[0];
	$mutation_hash{$line_array[0]}=1;
}

open INPUT, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_isolates.txt";
while(<INPUT>){
	$line=$_;
	chomp $line;
	
	push @isolate_array, $line;
}

open INPUT, "${base_dir}/Analysis/Mutation_intergenic_annotation/${species}_Mutation_intergenic_annotation/${species}_SNPs_intergenic.tab";
while(<INPUT>){
	$line=$_;
	chomp $line;
	@line_array=split(/\t/, $line, -1);
	
	if($line !~ /^Position\tCategory/){
		if($mutation_hash{$line_array[0]} && $line_array[1] eq "Promoter"){
			$mutation_ref_hash{$line_array[0]}=",$line_array[6],";
			$mutation_snp_hash{$line_array[0]}=",$line_array[7],";
		}
	}
}

foreach $mutation(@mutation_array){
	print OUTPUT "\t$mutation";
}
print OUTPUT "\n";

foreach $isolate(@isolate_array){
	print OUTPUT "$isolate";
	foreach $mutation(@mutation_array){
		if($mutation_ref_hash{$mutation} =~ /,$isolate,/){
			print OUTPUT "\tRef";
		}elsif($mutation_snp_hash{$mutation} =~ /,$isolate,/){
			print OUTPUT "\tSNP";
		}else{
			print OUTPUT "\tN";
		}
	}
	
	print OUTPUT "\n";
}

print "$species promoter heatmap created.\n";

print LOG "$species promoter heatmap created.\n";

