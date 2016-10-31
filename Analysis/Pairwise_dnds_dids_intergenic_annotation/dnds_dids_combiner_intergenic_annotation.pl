#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_dnds_dids_intergenic_annotation.csv";

open INPUT_DNDS, "${base_dir}/Analysis/Pairwise_dnds_dids/${species}_Pairwise_dnds_dids/${species}_dnds_half_a.txt";
while(<INPUT_DNDS>){
	$line=$_;
	chomp $line;
	if($line =~ /^(\S+)\s+(\S+)\s+(\S+.*)/){
		$isolate_1=$1;$isolate_2=$2;$info=$3;
		
		$info=~tr/\-/ /;
		$info=~tr/\(/ /;
		$info=~tr/\)/ /;
		
		if($info =~ /(\S+)\s+(\S+)\s+(\S+)/){
			$dnds=$1;$dn=$2;$ds=$3;
			
			push @isolate_1_array, $isolate_1;
			push @isolate_2_array, $isolate_2;
			push @dnds_array, $dnds;
			push @ds_div_array, $ds;
		}
	}
}

open INPUT_DNDS, "${base_dir}/Analysis/Pairwise_dnds_dids/${species}_Pairwise_dnds_dids/${species}_dnds_half_b.txt";
while(<INPUT_DNDS>){
	$line=$_;
	chomp $line;
	if($line =~ /(\S+)\s+(\S+)\s+(\S+.*)/){
		$isolate_1=$1;$isolate_2=$2;$info=$3;
		
		$info=~tr/\-/ /;
		$info=~tr/\(/ /;
		$info=~tr/\)/ /;
		
		if($info =~ /(\S+)\s+(\S+)\s+(\S+)/){
			$dnds=$1;$dn=$2;$ds=$3;
			
			push @ds_array, $ds;
		}
	}
}

open INPUT_DI, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_rbs_di.csv";
while(<INPUT_DI>){
	$line=$_;
	chomp $line;
	@line_array=split(/,/, $line, -1);
	if($line !~ /^Isolate/){
		push @di_rbs_array, $line_array[4];
	}
}

open INPUT_DI, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_promoter_di.csv";
while(<INPUT_DI>){
	$line=$_;
	chomp $line;
	@line_array=split(/,/, $line, -1);
	if($line !~ /^Isolate/){
		push @di_pro_array, $line_array[4];
	}
}

open INPUT_DI, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_terminator_di.csv";
while(<INPUT_DI>){
	$line=$_;
	chomp $line;
	@line_array=split(/,/, $line, -1);
	if($line !~ /^Isolate/){
		push @di_ter_array, $line_array[4];
	}
}

open INPUT_DI, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_non_coding_RNA_di.csv";
while(<INPUT_DI>){
	$line=$_;
	chomp $line;
	@line_array=split(/,/, $line, -1);
	if($line !~ /^Isolate/){
		push @di_rna_array, $line_array[4];
	}
}

open INPUT_DI, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_unannotated_di.csv";
while(<INPUT_DI>){
	$line=$_;
	chomp $line;
	@line_array=split(/,/, $line, -1);
	if($line !~ /^Isolate/){
		push @di_una_array, $line_array[4];
	}
}

$array_len=scalar(@di_una_array);

for($i=0; $i<$array_len; $i++){
	if($ds_div_array[$i] > 0 && $ds_array[$i] > 0){
		$dids_rbs_array[$i]=($di_rbs_array[$i]/$ds_div_array[$i]);
		$dids_pro_array[$i]=($di_pro_array[$i]/$ds_div_array[$i]);
		$dids_ter_array[$i]=($di_ter_array[$i]/$ds_div_array[$i]);
		$dids_rna_array[$i]=($di_rna_array[$i]/$ds_div_array[$i]);
		$dids_una_array[$i]=($di_una_array[$i]/$ds_div_array[$i]);
	}else{
		$dids_rbs_array[$i]="NA";
		$dids_pro_array[$i]="NA";
		$dids_ter_array[$i]="NA";
		$dids_rna_array[$i]="NA";
		$dids_una_array[$i]="NA";
		$zero_ds_hash{$i}=1;
	}
}

print OUTPUT "Isolate_1,Isolate_2,dN/dS,dI/dS_rbs,dI/dS_promoter,dI/dS_terminator,dI/dS_non_coding_RNA,dI/dS_unannotated,dS\n";

for($i=0; $i<$array_len; $i++){
	if(!$zero_ds_hash{$i}){
		print OUTPUT "$isolate_1_array[$i],$isolate_2_array[$i],$dnds_array[$i],$dids_rbs_array[$i],$dids_pro_array[$i],$dids_ter_array[$i],$dids_rna_array[$i],$dids_una_array[$i],$ds_array[$i]\n";
	}
}

print "$species intergenic annotation dids combined.\n";

print LOG "$species intergenic annotation dids combined.\n";

