#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

open OUTPUT_G, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_gene_coordinates.tab";
print OUTPUT_G "Name\tStart\tEnd\tLength\tType\tStrand\n";
open OUTPUT_I, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_intergenic_coordinates.tab";
print OUTPUT_I "Name\tStart\tEnd\tLength\tType\n";

open INPUT, "${base_dir}/Data/GFF_files/$species.gff";
while(<INPUT>){
	$line=$_;
	chomp $line;
	@line_array=split(/\t/, $line, -1);
	
	$sta="";
	$end="";
	$strand="";
	$id="";
	if($line !~ /^##/){
		if($line_array[2] ne "sig_peptide" && $line_array[2] ne "misc_RNA"){
			$sta=$line_array[3];
			$end=$line_array[4];
			$type=$line_array[2];
		
			if($line_array[6] eq "+"){
				$strand="Forward";
			}elsif($line_array[6] eq "-"){
				$strand="Reverse";
			}
		
			$len=(($end - $sta) + 1);
		
			if($line_array[8] =~ /ID=(${species}_[^;]+);/){
				$id=$1;
			}
		
			@tmp_array=();
			@tmp_array=("$id", "$sta", "$end", "$len", "$type", "$strand");
		
			push @gene_array, [@tmp_array];
		
			print OUTPUT_G "$id\t$sta\t$end\t$len\t$type\t$strand\n";
		
		}
	}elsif($line =~ /^##sequence-region\s+\S+\s+(\d+)\s+(\d+)/){
		$seq_sta=$1;
		$seq_end=$2;
	}elsif($line =~ /^##FASTA/){
		last;
	}
}

$gene_count=scalar(@gene_array);

$int_count=0;
for($i=0; $i<$gene_count; $i++){
	if($i == 0){
		$int_sta=$seq_sta;
		$int_end=($gene_array[$i][1] - 1);
		
		if($gene_array[$i][5] eq "Forward"){
			$int_type="Intergenic_co-oriented_F";
		}elsif($gene_array[$i][5] eq "Reverse"){
			$int_type="Intergenic_double_terminator";
		}
	}elsif($i == ($gene_count - 1)){
		$int_sta=($gene_array[($i-1)][2] + 1);
		$int_end=$seq_end;
		
		if($gene_array[($i-1)][5] eq "Forward"){
			$int_type="Intergenic_co-oriented_F";
		}elsif($gene_array[($i-1)][5] eq "Reverse"){
			$int_type="Intergenic_double_terminator";
		}
	}else{
		$int_sta=($gene_array[($i-1)][2] + 1);
		$int_end=($gene_array[$i][1] - 1);
		
		if($gene_array[($i-1)][5] eq "Forward" && $gene_array[$i][5] eq "Forward"){
			$int_type="Intergenic_co-oriented_F";
		}elsif($gene_array[($i-1)][5] eq "Forward" && $gene_array[$i][5] eq "Reverse"){
			$int_type="Intergenic_double_terminator";
		}elsif($gene_array[($i-1)][5] eq "Reverse" && $gene_array[$i][5] eq "Forward"){
			$int_type="Intergenic_double_promoter";
		}elsif($gene_array[($i-1)][5] eq "Reverse" && $gene_array[$i][5] eq "Reverse"){
			$int_type="Intergenic_co-oriented_R";
		}
	}
	
	$int_len=(($int_end - $int_sta) + 1);
	
	if($int_len > 0){
		$int_count++;
		
		$int_id="${species}_intergenic_$int_count";
		
		print OUTPUT_I "$int_id\t$int_sta\t$int_end\t$int_len\t$int_type\n";
	}
}

print "$species gene intergenic coordinates extracted.\n";

print LOG "$species gene intergenic coordinates extracted.\n";

