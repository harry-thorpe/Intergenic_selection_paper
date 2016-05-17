#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_site_counts.csv";

print OUTPUT "Category,A,T,G,C,Total,GC_content\n";

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

$comp_hash{"A"}="T";
$comp_hash{"T"}="A";
$comp_hash{"G"}="C";
$comp_hash{"C"}="G";

@base_array=("A", "T", "G", "C");

open INPUT, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_core_gene_alignment.fasta";
while(<INPUT>){
	$line=$_;
	chomp $line;
	
	if($line =~ /^([ATGCN]+)/){
		$seq=$1;
		$seq_len=length($seq);
		@seq_array=split(//, $seq);
		
		for($i=0; $i<$seq_len; $i+=3){
		
			$refcod="$seq_array[$i]$seq_array[($i+1)]$seq_array[($i+2)]";
		
			foreach $base(@base_array){
				if($base ne $seq_array[$i]){
					$alicod="$base$seq_array[($i+1)]$seq_array[($i+2)]";
				
					if($cod{$refcod} eq $cod{$alicod}){
						$site_hash{"Synonymous"}{$seq_array[$i]}++;
					}elsif($cod{$refcod} ne $cod{$alicod}){
						if($cod{$alicod} ne "STO"){
							$site_hash{"Nonsynonymous"}{$seq_array[$i]}++;
						}elsif($cod{$alicod} eq "STO"){
							$site_hash{"Nonsense"}{$seq_array[$i]}++;
						}
					}
				}
			}
		
			foreach $base(@base_array){
				if($base ne $seq_array[($i+1)]){
					$alicod="$seq_array[$i]$base$seq_array[($i+2)]";
				
					if($cod{$refcod} eq $cod{$alicod}){
						$site_hash{"Synonymous"}{$seq_array[($i+1)]}++;
					}elsif($cod{$refcod} ne $cod{$alicod}){
						if($cod{$alicod} ne "STO"){
							$site_hash{"Nonsynonymous"}{$seq_array[($i+1)]}++;
						}elsif($cod{$alicod} eq "STO"){
							$site_hash{"Nonsense"}{$seq_array[($i+1)]}++;
						}
					}
				}
			}
		
			foreach $base(@base_array){
				if($base ne $seq_array[($i+2)]){
					$alicod="$seq_array[$i]$seq_array[($i+1)]$base";
				
					if($cod{$refcod} eq $cod{$alicod}){
						$site_hash{"Synonymous"}{$seq_array[($i+2)]}++;
					}elsif($cod{$refcod} ne $cod{$alicod}){
						if($cod{$alicod} ne "STO"){
							$site_hash{"Nonsynonymous"}{$seq_array[($i+2)]}++;
						}elsif($cod{$alicod} eq "STO"){
							$site_hash{"Nonsense"}{$seq_array[($i+2)]}++;
						}
					}
				}
			}
		}
		
		last;
	}
}

open INPUT, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_core_intergenic_alignment.fasta";
while(<INPUT>){
	$line=$_;
	chomp $line;
	
	if($line =~ /^([ATGCN]+)/){
		$seq=$1;
		$seq_len=length($seq);
		@seq_array=split(//, $seq);
		
		for($i=0; $i<$seq_len; $i++){
			$site_hash{"Intergenic"}{$seq_array[$i]}++;
		}
		
		last;
	}
}

@type_array=("Synonymous", "Nonsynonymous", "Nonsense", "Intergenic");

foreach $type(@type_array){
	$a=$site_hash{$type}{"A"};
	$t=$site_hash{$type}{"T"};
	$g=$site_hash{$type}{"G"};
	$c=$site_hash{$type}{"C"};
	
	if($type ne "Intergenic"){
		$a=($a/3);
		$t=($t/3);
		$g=($g/3);
		$c=($c/3);
	}
	$tot=($a + $t + $g + $c);
	$gc_con=(($g + $c) / $tot);
	
	print OUTPUT "$type,$a,$t,$g,$c,$tot,$gc_con\n";
}


print "$species sites counted.\n";

print LOG "$species sites counted.\n";

