#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

@category_array=("terminator_stem", "terminator_loop");

$category_count=scalar(@category_array);

for($category_i=0; $category_i<$category_count; $category_i++){
	
	$category=$category_array[$category_i];
	
	mkdir "${base_dir}/Analysis/${analysis}/${species}_${analysis}/Intergenic_${category}_files";

	$count=0;
	open INPUT, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_core_intergenic_${category}_alignment.fasta";
	while(<INPUT>){
		$line=$_;
		chomp $line;
	
		if($line =~/^>(\S+)/){
			$id=$1;
		}elsif($line =~ /^([ATGCN]+)/){
			$seq=$1;
		
			open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/Intergenic_${category}_files/$id.fasta";
			print OUTPUT ">$id\n$seq\n";
		
			$count++;
			#print "Isolate $count completed.\n";
		}
	}
}

print "$species terminator stem loop alignments split.\n";

print LOG "$species terminator stem loop alignments split.\n";

