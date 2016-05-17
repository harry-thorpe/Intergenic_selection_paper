#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_core_gene_alignment_no_STO.fasta";

%remove_sites_hash=();

$count=0;
open INGEN, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_core_gene_alignment.fasta";
while(<INGEN>){
	if(/^>(\S+)/){
		$id=$1;
	}elsif(/^([ATGCN]+)/){
		$seq=$1;
		@seq_array=split(//, $seq);
		$len=scalar(@seq_array);
		for($i=0; $i<$len; $i+=3){
			$codon=$seq_array[$i].$seq_array[($i+1)].$seq_array[($i+2)];
			if($codon=~/(TAG|TGA|TAA)/){
				$remove_sites_hash{$i}++;
			}
		}
		$count++;
		#print "$count isolates read.\n";
	}
}

$count=0;
open INGEN, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_core_gene_alignment.fasta";
while(<INGEN>){
	if(/^>(\S+)/){
		$id=$1;
		print OUTPUT ">$id\n";
	}elsif(/^([ATGCN]+)/){
		$seq=$1;
		@seq_array=split(//, $seq);
		$len=scalar(@seq_array);
		for($i=0; $i<$len; $i+=3){
			$codon=$seq_array[$i].$seq_array[($i+1)].$seq_array[($i+2)];
			if(!$remove_sites_hash{$i}){
				print OUTPUT "$codon";
			}
		}
		print OUTPUT "\n";
		$count++;
		#print "$count isolates completed.\n";
	}
}

print "$species core gene stop codons removed.\n";

print LOG "$species core gene stop codons removed.\n";

