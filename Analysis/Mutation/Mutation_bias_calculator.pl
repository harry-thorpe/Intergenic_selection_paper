#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_mutation_bias.tab";

print OUTPUT "Category\tMutation\tCount\n";

$comp_hash{"A"}="T";
$comp_hash{"T"}="A";
$comp_hash{"G"}="C";
$comp_hash{"C"}="G";

$cum_a=0;
$cum_t=0;
$cum_g=0;
$cum_c=0;
$cum_total=0;
open INPUT, "${base_dir}/Analysis/Sequence_summary/${species}_Sequence_summary/${species}_site_counts.csv";
while(<INPUT>){
	$line=$_;
	chomp $line;
	@line_array=split(/,/, $line);
	
	if($line !~ /^Category,A,T,G,C,Total,GC_content/){
		$cat=$line_array[0];
		$a=$line_array[1];
		$t=$line_array[2];
		$g=$line_array[3];
		$c=$line_array[4];
		$total=$line_array[5];
		
		$at=(($a + $t) / $total);
		$gc=(($g + $c) / $total);
		
		$gc_hash{$cat}{"AT"}=$at;
		$gc_hash{$cat}{"AG"}=$at;
		$gc_hash{$cat}{"AC"}=$at;
		
		$gc_hash{$cat}{"TA"}=$at;
		$gc_hash{$cat}{"TG"}=$at;
		$gc_hash{$cat}{"TC"}=$at;
		
		$gc_hash{$cat}{"CA"}=$gc;
		$gc_hash{$cat}{"CT"}=$gc;
		$gc_hash{$cat}{"CG"}=$gc;
		
		$gc_hash{$cat}{"GA"}=$gc;
		$gc_hash{$cat}{"GT"}=$gc;
		$gc_hash{$cat}{"GC"}=$gc;
		
		$cum_a=($cum_a+$a);
		$cum_t=($cum_t+$t);
		$cum_g=($cum_g+$g);
		$cum_c=($cum_c+$c);
		$cum_total=($cum_total+$total);
	}
}

$cum_at=(($cum_a + $cum_t) / $cum_total);
$cum_gc=(($cum_g + $cum_c) / $cum_total);

$gc_hash{"All"}{"AT"}=$cum_at;
$gc_hash{"All"}{"AG"}=$cum_at;
$gc_hash{"All"}{"AC"}=$cum_at;

$gc_hash{"All"}{"TA"}=$cum_at;
$gc_hash{"All"}{"TG"}=$cum_at;
$gc_hash{"All"}{"TC"}=$cum_at;

$gc_hash{"All"}{"CA"}=$cum_gc;
$gc_hash{"All"}{"CT"}=$cum_gc;
$gc_hash{"All"}{"CG"}=$cum_gc;

$gc_hash{"All"}{"GA"}=$cum_gc;
$gc_hash{"All"}{"GT"}=$cum_gc;
$gc_hash{"All"}{"GC"}=$cum_gc;

open INPUT, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_mutations.tab";
while(<INPUT>){
	$line=$_;
	chomp $line;
	@line_array=split(/\t/, $line);
	
	if($line !~ /^Position\tCategory/){
		$cat=$line_array[1];
		$ref=$line_array[2];
		$snp=$line_array[3];
		$freq=$line_array[5];
		
		if($freq == 1){
			@ref_array=split(//, $ref);
			@snp_array=split(//, $snp);
			
			$len=scalar(@ref_array);
			
			$ref_base="";
			$snp_base="";
			$snp_count=0;
			for($i=0; $i<$len; $i++){
				if($ref_array[$i] ne $snp_array[$i]){
					$snp_count++;
					
					$ref_base=$ref_array[$i];
					$snp_base=$snp_array[$i];
				}
			}
			
			if($snp_count == 1){
				$snp_type="$ref_base$snp_base";
				$snp_type_comp="$comp_hash{$ref_base}$comp_hash{$snp_base}";
				
				@snp_type_array=("$snp_type", "$snp_type_comp");
				@snp_type_array=sort(@snp_type_array);
				
				$snp_hash{$cat}{$snp_type_array[0]}++;
				$snp_hash{"All"}{$snp_type_array[0]}++;
				$snp_hash{$cat}{$snp_type_array[1]}++;
				$snp_hash{"All"}{$snp_type_array[1]}++;
				#$snp_count_hash{$cat}++;
				#$snp_count_hash{"All"}++;
				
			}else{
				print "bad SNP.\n";
			}
		}
	}
}

@cat_array=keys(%snp_hash);
@cat_array=sort(@cat_array);

@type_array=keys(%{$snp_hash{$cat_array[0]}});
@type_array=sort(@type_array);

foreach $cat(@cat_array){
	foreach $type(@type_array){
		if($snp_hash{$cat}{$type}){
			$prop=($snp_hash{$cat}{$type}/$gc_hash{$cat}{$type});
		}else{
			$prop=0;
		}
		
		$snp_gc_corrected_hash{$cat}{$type}=$prop;
		$snp_gc_corrected_total_hash{$cat}+=$prop;
	}
}

foreach $cat(@cat_array){
	foreach $type(@type_array){
		$rel_prop=int(($snp_gc_corrected_hash{$cat}{$type}/$snp_gc_corrected_total_hash{$cat}) * 100);
		
		print OUTPUT "$cat\t$type\t$rel_prop\n";
	}
}

print "$species mutation bias calculated.\n";

print LOG "$species mutation bias calculated.\n";

