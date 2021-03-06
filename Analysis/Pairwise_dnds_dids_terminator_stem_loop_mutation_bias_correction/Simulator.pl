#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

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

open LOG, ">>${base_dir}/Analysis/log.txt";

$in_gene_file="${base_dir}/Analysis/Pairwise_dnds_dids/${species}_Pairwise_dnds_dids/${species}_core_gene_alignment_no_STO.fasta";
$in_int_ter_stem_file="${base_dir}/Analysis/Pairwise_dnds_dids_terminator_stem_loop/${species}_Pairwise_dnds_dids_terminator_stem_loop/${species}_core_intergenic_terminator_stem_alignment.fasta";
$in_int_ter_loop_file="${base_dir}/Analysis/Pairwise_dnds_dids_terminator_stem_loop/${species}_Pairwise_dnds_dids_terminator_stem_loop/${species}_core_intergenic_terminator_loop_alignment.fasta";

open OUTPUT_GENE, ">${base_dir}/Analysis/${analysis}/${species}_$analysis/${species}_core_gene_alignment_simulated.fasta";
open OUTPUT_INT_TER_STEM, ">${base_dir}/Analysis/${analysis}/${species}_$analysis/${species}_core_intergenic_terminator_stem_alignment_simulated.fasta";
open OUTPUT_INT_TER_LOOP, ">${base_dir}/Analysis/${analysis}/${species}_$analysis/${species}_core_intergenic_terminator_loop_alignment_simulated.fasta";

#$inmut="${base_dir}/Analysis/${analysis}/Mutation_bias.tab";
#$inmut="${base_dir}/Analysis/${analysis}/Mutation_no_bias.tab";
$inmut="${base_dir}/Analysis/Mutation/${species}_Mutation/${species}_mutation_bias.tab";

$strs=10;

$typecount=0;
open INMUT, $inmut;
while(<INMUT>){
	if(/^All\s+(\S+)\s+(\S+)/){
		$type=$1;$count=$2;
		for $x($typecount..($typecount+($count-1))){
			$mut{$x}=$type;
		}
		$typecount=($typecount+$count);
	}
}
$tmuts=$typecount;

open INPUT, $in_gene_file;
while(<INPUT>){
	if(/^([ATGCN]+)/){
		$gene_seq=$1;
		$gene_len=length($gene_seq);
		last;
	}
}
open INPUT, $in_int_ter_stem_file;
while(<INPUT>){
	if(/^([ATGCN]+)/){
		$int_ter_stem_seq=$1;
		$int_ter_stem_len=length($int_ter_stem_seq);
		last;
	}
}
open INPUT, $in_int_ter_loop_file;
while(<INPUT>){
	if(/^([ATGCN]+)/){
		$int_ter_loop_seq=$1;
		$int_ter_loop_len=length($int_ter_loop_seq);
		last;
	}
}

$seq="$gene_seq$int_ter_stem_seq$int_ter_loop_seq";
$tsites=length($seq);

$max_muts=int(($tsites / 100) / 2);

@max_muts_array=(int($max_muts*0.2), int($max_muts*0.4), int($max_muts*0.6), int($max_muts*0.8), $max_muts);

#print "@max_muts_array\n";

foreach $max_muts(@max_muts_array){
	for $y(1..$strs){
	
		$muts=0;
		@mutseq=split(//, $seq);
	
		while($muts < $max_muts){
		
			$ransite=int(rand($tsites));
			$ranmut=int(rand($tmuts));
			@mutbase=split(//, $mut{$ranmut});
		
			if($mutbase[0] eq $mutseq[$ransite]){
				if($ransite < $gene_len){
					if($ransite%3==0){
						$a=$ransite;
						$refcod="$mutseq[$a]$mutseq[($a+1)]$mutseq[($a+2)]";
						$mutcod="$mutbase[1]$mutseq[($a+1)]$mutseq[($a+2)]";
					}elsif($ransite%3==1){
						$a=($ransite-1);
						$refcod="$mutseq[$a]$mutseq[($a+1)]$mutseq[($a+2)]";
						$mutcod="$mutseq[$a]$mutbase[1]$mutseq[($a+2)]";
					}elsif($ransite%3==2){
						$a=($ransite-2);
						$refcod="$mutseq[$a]$mutseq[($a+1)]$mutseq[($a+2)]";
						$mutcod="$mutseq[$a]$mutseq[($a+1)]$mutbase[1]";
					}
					if($cod{$refcod} eq "STO" or $cod{$mutcod} eq "STO"){
						#print "STO $refcod to $mutcod\n";
					}elsif($cod{$refcod} ne "STO" && $cod{$mutcod} ne "STO"){
						$mutseq[$ransite]=$mutbase[1];
						$muts++;
						#print "Strain $y, mutation $muts completed.\n";
					}
				}else{
					$mutseq[$ransite]=$mutbase[1];
					$muts++;
					#print "Strain $y, mutation $muts completed.\n";
				}
			}
		}
		
		#print "Strain $y, mutation $muts completed.\n";
		
		$min=0;
		$max=$gene_len;
		print OUTPUT_GENE ">${y}_$max_muts\n";
		for ($i=$min; $i<$max; $i++){
			print OUTPUT_GENE "$mutseq[$i]";
		}
		print OUTPUT_GENE "\n";
		
		$min=$max;
		$max=$max+$int_ter_stem_len;
		print OUTPUT_INT_TER_STEM ">${y}_$max_muts\n";
		for ($i=$min; $i<$max; $i++){
			print OUTPUT_INT_TER_STEM "$mutseq[$i]";
		}
		print OUTPUT_INT_TER_STEM "\n";
		
		$min=$max;
		$max=$max+$int_ter_loop_len;
		print OUTPUT_INT_TER_LOOP ">${y}_$max_muts\n";

		for ($i=$min; $i<$max; $i++){
			print OUTPUT_INT_TER_LOOP "$mutseq[$i]";
		}
		print OUTPUT_INT_TER_LOOP "\n";
	}
}

print "$species dnds dids simulated.\n";

print LOG "$species dnds dids simulated.\n";

