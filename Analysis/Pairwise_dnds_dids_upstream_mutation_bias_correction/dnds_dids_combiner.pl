#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_dnds_dids_intergenic_upstream_simulated.csv";

open INPUT_DNDS, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_dnds_simulated.csv";
while(<INPUT_DNDS>){
	$line=$_;
	chomp $line;
	@line_array=split(/,/, $line);
	
	if($line !~ /^Isolate_1,Isolate_2/){
		push @isolate_1_array, $line_array[0];
		push @isolate_2_array, $line_array[1];
		#push @dn_array, $line_array[2];
		push @ds_array, $line_array[3];
		push @dnds_array, $line_array[4];
	}
}

open INPUT_DI, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_all_di_simulated.csv";
while(<INPUT_DI>){
	$line=$_;
	chomp $line;
	@line_array=split(/,/, $line);
	
	if($line !~ /^Isolate_1,Isolate_2/){
		push @di_all_array, $line_array[4];
	}
}

open INPUT_DI, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_upstream_30_di_simulated.csv";
while(<INPUT_DI>){
	$line=$_;
	chomp $line;
	@line_array=split(/,/, $line);
	
	if($line !~ /^Isolate_1,Isolate_2/){
		push @di_ups_array, $line_array[4];
	}
}

$array_len=scalar(@di_all_array);

for($i=0; $i<$array_len; $i++){
	if($ds_array[$i] > 0){
		$dids_all_array[$i]=($di_all_array[$i]/$ds_array[$i]);
		$dids_ups_array[$i]=($di_ups_array[$i]/$ds_array[$i]);
	}else{
		$dids_all_array[$i]="NA";
		$dids_ups_array[$i]="NA";
		$zero_ds_hash{$i}=1;
	}
}

print OUTPUT "Isolate_1,Isolate_2,dN/dS,dI/dS,dI/dS_upstream_30,dS\n";

for($i=0; $i<$array_len; $i++){
	if(!$zero_ds_hash{$i}){
		print OUTPUT "$isolate_1_array[$i],$isolate_2_array[$i],$dnds_array[$i],$dids_all_array[$i],$dids_ups_array[$i],$ds_array[$i]\n";
	}
}

print "$species dnds dids results combined.\n";

print LOG "$species dnds dids results combined.\n";

