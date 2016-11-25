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

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}$threshold_folder/${species}_mutations.tab";

print OUTPUT "Position\tCategory\tReference\tSNP\tReference_count\tSNP_count\tReference_isolates\tSNP_isolates\n";

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

$count=0;
open INPUT, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment$threshold_folder/${species}_core_gene_alignment.fasta";
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
				
				$codon_snp_pos=$i;
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
			
			$biallelic_codon_hash{$codon_pos}=$mut;
		}
	}
}

%codon_hash=();

@biallelic_codon_array=keys(%biallelic_codon_hash);
@biallelic_codon_array=sort { $a <=> $b } @biallelic_codon_array;

$count=0;
open INPUT, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment$threshold_folder/${species}_core_gene_alignment.fasta";
while(<INPUT>){
	$line=$_;
	chomp $line;
	
	if($line =~ /^>(\S+)/){
		$id=$1;
		
	}elsif($line =~ /^([ATGCN]+)/){
		@seq_array=split(//, $line);
		
		foreach $codon_pos(@biallelic_codon_array){

			$codon="$seq_array[$codon_pos]$seq_array[($codon_pos+1)]$seq_array[($codon_pos+2)]";
			
			if($codon !~ /N/){
				$codon_hash{$codon_pos}{$codon}++;
				
				if(!$codon_isolate_hash{$codon_pos}{$codon}){
					$codon_isolate_hash{$codon_pos}{$codon}=$id;
				}else{
					$codon_isolate_hash{$codon_pos}{$codon}="$codon_isolate_hash{$codon_pos}{$codon},$id";
				}
			}
		}
		
		$count++;
		#print "Isolate $count completed.\n";
	}
}

foreach $codon_pos(@biallelic_codon_array){
	@codon_array=sort { $codon_hash{$codon_pos}{$b} <=> $codon_hash{$codon_pos}{$a} } keys(%{$codon_hash{$codon_pos}});
	
	$codon_count=scalar(@codon_array);
	
	if($codon_count == 2){
		
		@codon_1_array=split(//, $codon_array[0]);
		@codon_2_array=split(//, $codon_array[1]);
		
		$codon_snp=0;
		for($i=0; $i<3; $i++){
			if($codon_1_array[$i] ne $codon_2_array[$i]){
				$codon_snp++;
				
				$codon_snp_pos=$i;
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
			
			$real_base_pos=($codon_pos + $codon_snp_pos + 1);
			
			print OUTPUT "$real_base_pos\t$mut\t$codon_array[0]\t$codon_array[1]\t$codon_hash{$codon_pos}{$codon_array[0]}\t$codon_hash{$codon_pos}{$codon_array[1]}\t$codon_isolate_hash{$codon_pos}{$codon_array[0]}\t$codon_isolate_hash{$codon_pos}{$codon_array[1]}\n";
		}
	}
}


$count=0;
open INPUT, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment$threshold_folder/${species}_core_intergenic_alignment.fasta";
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
		
		$biallelic_base_hash{$base_pos}=$mut;
	}
}

%base_hash=();

@biallelic_base_array=keys(%biallelic_base_hash);
@biallelic_base_array=sort { $a <=> $b } @biallelic_base_array;

$count=0;
open INPUT, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment$threshold_folder/${species}_core_intergenic_alignment.fasta";
while(<INPUT>){
	$line=$_;
	chomp $line;
	
	if($line =~ /^>(\S+)/){
		$id=$1;
		
	}elsif($line =~ /^([ATGCN]+)/){
		@seq_array=split(//, $line);
		
		foreach $base_pos(@biallelic_base_array){
			if($seq_array[$base_pos] !~ /N/){
				$base_hash{$base_pos}{$seq_array[$base_pos]}++;
			
				if(!$base_isolate_hash{$base_pos}{$seq_array[$base_pos]}){
					$base_isolate_hash{$base_pos}{$seq_array[$base_pos]}=$id;
				}else{
					$base_isolate_hash{$base_pos}{$seq_array[$base_pos]}="$base_isolate_hash{$base_pos}{$seq_array[$base_pos]},$id";
				}
			}
		}
		
		$count++;
		#print "Isolate $count completed.\n";
	}
}

foreach $base_pos(@biallelic_base_array){
	@base_array=sort { $base_hash{$base_pos}{$b} <=> $base_hash{$base_pos}{$a} } keys(%{$base_hash{$base_pos}});
	
	$base_count=scalar(@base_array);
	
	if($base_count == 2){
		
		$mut="Intergenic";
		
		$real_base_pos=($base_pos + 1);
		
		print OUTPUT "$real_base_pos\t$mut\t$base_array[0]\t$base_array[1]\t$base_hash{$base_pos}{$base_array[0]}\t$base_hash{$base_pos}{$base_array[1]}\t$base_isolate_hash{$base_pos}{$base_array[0]}\t$base_isolate_hash{$base_pos}{$base_array[1]}\n";
	}
}

print "$species mutations called.\n";

print LOG "$species mutations called.\n";

