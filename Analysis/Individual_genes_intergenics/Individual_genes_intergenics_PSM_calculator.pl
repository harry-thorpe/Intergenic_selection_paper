#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_core_gene_intergenic_PSM.tab";

print OUTPUT "Region\tCategory\tSNPs\tSingletons\tPSM\n";

$cod{TTT}="Phe";  $cod{TCT}="Ser";  $cod{TAT}="Tyr";  $cod{TGT}="Cys";
$cod{TTC}="Phe";  $cod{TCC}="Ser";  $cod{TAC}="Tyr";  $cod{TGC}="Cys";
$cod{TTA}="Leu";  $cod{TCA}="Ser";  $cod{TAA}="STO";  $cod{TGA}="STO";
$cod{TTG}="Leu";  $cod{TCG}="Ser";  $cod{TAG}="STO";  $cod{TGG}="Trp";

$cod{CTT}="Leu";  $cod{CCT}="Pro";  $cod{CAT}="His";  $cod{CGT}="Arg";
$cod{CTC}="Leu";  $cod{CCC}="Pro";  $cod{CAC}="His";  $cod{CGC}="Arg";
$cod{CTA}="Leu";  $cod{CCA}="Pro";  $cod{CAA}="Gln";  $cod{CGA}="Arg";
$cod{CTG}="Leu";  $cod{CCG}="Pro";  $cod{CAG}="Gln";  $cod{CGG}="Arg";

$cod{ATT}="Ile";  $cod{ACT}="Thr";  $cod{AAT}="Asn";  $cod{AGT}="Ser";
$cod{ATC}="Ile";  $cod{ACC}="Thr";  $cod{AAC}="Asn";  $cod{AGC}="Ser";
$cod{ATA}="Ile";  $cod{ACA}="Thr";  $cod{AAA}="Lys";  $cod{AGA}="Arg";
$cod{ATG}="Met";  $cod{ACG}="Thr";  $cod{AAG}="Lys";  $cod{AGG}="Arg";

$cod{GTT}="Val";  $cod{GCT}="Ala";  $cod{GAT}="Asp";  $cod{GGT}="Gly";
$cod{GTC}="Val";  $cod{GCC}="Ala";  $cod{GAC}="Asp";  $cod{GGC}="Gly";
$cod{GTA}="Val";  $cod{GCA}="Ala";  $cod{GAA}="Glu";  $cod{GGA}="Gly";
$cod{GTG}="Val";  $cod{GCG}="Ala";  $cod{GAG}="Glu";  $cod{GGG}="Gly";

open INCOOR, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_core_genes.tab";
while(<INCOOR>){
	$line=$_;
	chomp $line;
	if($line =~ /^(\S+)\s+\d+\s+\d+\s+\d+\s+\S+\s+\S+/){
		$gene_id=$1;
		if($line !~ /^Name\tStart\tEnd\tLength\tType/){
			push @gene_array, $gene_id;
		}
	}
}

@mut_array=("Nonsynonymous", "Synonymous", "Nonsense");

foreach $gene(@gene_array){
	%codon_hash=();
	%mut_type_hash=();
	%mut_type_singleton_hash=();
	
	$count=0;
	open INPUT, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_gene_files/${gene}.fasta";
	while(<INPUT>){
		$line=$_;
		chomp $line;
	
		if($line =~ /^>(\S+)/){
			$id=$1;
		
		}elsif($line =~ /^([ATGCN]+)/){
			@seq_array=split(//, $line);
		
			$seq_len=length($line);
		
			for($i=0; $i<$seq_len; $i+=3){
				$codon="$seq_array[$i]$seq_array[($i+1)]$seq_array[($i+2)]";
			
				if($codon !~ /N/){
					$codon_hash{$i}{$codon}++;
				}
			}
		
			$count++;
			#print "Isolate $count completed.\n";
		}
	}

	@codon_pos_array=keys(%codon_hash);
	@codon_pos_array=sort { $a <=> $b } @codon_pos_array;

	foreach $codon_pos(@codon_pos_array){
		@codon_array=sort { $codon_hash{$codon_pos}{$b} <=> $codon_hash{$codon_pos}{$a} } keys(%{$codon_hash{$codon_pos}});
	
		$codon_count=scalar(@codon_array);
	
		if($codon_count == 2){
		
			@codon_1_array=split(//, $codon_array[0]);
			@codon_2_array=split(//, $codon_array[1]);
		
			$codon_snp=0;
			for($i=0; $i<3; $i++){
				if($codon_1_array[$i] ne $codon_2_array[$i]){
					$codon_snp++;
				
					#$codon_snp_pos=$i;
				}
			}
		
			if($codon_snp == 1){
				if($cod{$codon_array[0]} eq $cod{$codon_array[1]}){
					$mut="Synonymous";
				}elsif($cod{$codon_array[0]} ne $cod{$codon_array[1]}){
					if($cod{$codon_array[1]} =~ /STO/){
						$mut="Nonsense";
					}else{
						$mut="Nonsynonymous";
					}
				}
				
				$mut_type_hash{$mut}++;
				
				if($codon_hash{$codon_pos}{$codon_array[1]} == 1){
					$mut_type_singleton_hash{$mut}++;
				}
			}
		}
	}
	
	foreach $mut(@mut_array){
		if($mut_type_hash{$mut} && $mut_type_singleton_hash{$mut}){
			if($mut_type_hash{$mut} > 0){
				$psm=($mut_type_singleton_hash{$mut} / $mut_type_hash{$mut});
	
				print OUTPUT "$gene\t$mut\t$mut_type_hash{$mut}\t$mut_type_singleton_hash{$mut}\t$psm\n";
			}
		}
	}
}

open IN_COOR, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_core_intergenics.tab";
while(<IN_COOR>){
	$line=$_;
	chomp $line;
	if($line =~ /^(\S+)\s+\d+\s+\d+\s+\d+\s+\S+/){
		$int_id=$1;
		if($line !~ /^Name\tStart\tEnd\tLength\tType/){
				push @intergenic_array, $int_id;
		}
	}

}

foreach $intergenic(@intergenic_array){
	%base_hash=();
	%mut_type_hash=();
	%mut_type_singleton_hash=();
	
	$count=0;
	open INPUT, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_intergenic_files/${intergenic}.fasta";
	while(<INPUT>){
		$line=$_;
		chomp $line;
	
		if($line =~ /^>(\S+)/){
			$id=$1;
		
		}elsif($line =~ /^([ATGCN]+)/){
	
			@seq_array=split(//, $line);
			
			$seq_len=length($line);
	
			for($i=0; $i<$seq_len; $i++){
				if($seq_array[$i] !~ /N/){
			
					$base_hash{$i}{$seq_array[$i]}++;
				}
			}
	
			$count++;
			#print "Isolate $count completed.\n";
		}
	}
	
	@base_pos_array=keys(%base_hash);
	@base_pos_array=sort { $a <=> $b } @base_pos_array;


	foreach $base_pos(@base_pos_array){
		@base_array=sort { $base_hash{$base_pos}{$b} <=> $base_hash{$base_pos}{$a} } keys(%{$base_hash{$base_pos}});

		$base_count=scalar(@base_array);

		if($base_count == 2){
	
			$mut="Intergenic";
	
			$mut_type_hash{$mut}++;
	
			if($base_hash{$base_pos}{$base_array[1]} == 1){
				$mut_type_singleton_hash{$mut}++;
			}
		}
	}
	
	$mut="Intergenic";
	
	if($mut_type_hash{$mut} && $mut_type_singleton_hash{$mut}){
		if($mut_type_hash{$mut} > 0){
			$psm=($mut_type_singleton_hash{$mut} / $mut_type_hash{$mut});
	
			print OUTPUT "$intergenic\t$mut\t$mut_type_hash{$mut}\t$mut_type_singleton_hash{$mut}\t$psm\n";
		}
	}
}

print "$species individual genes intergenics PSM calculated.\n";

print LOG "$species individual genes intergenics PSM calculated.\n";

