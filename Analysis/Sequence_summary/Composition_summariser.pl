#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_composition_summary.csv";

print OUTPUT "Triplet,Gene,Intergenic\n";

open OUTPUT2, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_GC_content_summary.csv";

print OUTPUT2 "Category,GC_content\n";

open OUTPUT3, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_reference_GC_content.txt";

%trp_gene_hash=();
$trp_gene=0;
open REF, "${base_dir}/Data/Reference_files/${species}_reference.fasta";
while(<REF>){
	if(/^([ATGCN]+)/){
		$seq=$1;
		@seq_array=split(//, $seq);
		$seq_len=length($seq);
		
		$n=$seq=~tr/N/N/;
		#print "$n\n";
		
		$gc=0;
		$atgc=0;
		for($i=0; $i<$seq_len; $i++){
			$base=$seq_array[$i];
	
			if($base =~ /G|C/){
				$gc++;
			}
			if($base =~ /A|T|G|C/){
				$atgc++;
			}
		}

		$gc_content=($gc / $atgc);
		
		print OUTPUT3 "$gc_content";
		
		open INPUT, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_core_genes.tab";
		while(<INPUT>){
			$line=$_;
			chomp $line;
			if($line =~ /^\S+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+\S+/){
				$beg=$1;$end=$2;$len=$3;$type=$4;
				if($line !~ /^Name\tStart\tEnd\tLength\tType/){
					$indbeg=($beg-1);
					$indend=($end-1);
				
					$type="Gene";
				
					for($i=$indbeg; $i<=($indend-2); $i++){
						$trp="$seq_array[$i]$seq_array[$i+1]$seq_array[$i+2]";
						$trp_rev=reverse($trp);
					
						@trp_array=($trp, $trp_rev);
						@trp_array=sort(@trp_array);
					
						$trp_gene_hash{$trp_array[0]}++;
						$trp_gene++;
					}
				
					$gc=0;
					$atgc=0;
					for($i=$indbeg; $i<=$indend; $i++){
						$base=$seq_array[$i];
				
						if($base =~ /G|C/){
							$gc++;
						}
						if($base =~ /A|T|G|C/){
							$atgc++;
						}
					}
			
					$gc_content=($gc / $atgc);
				
					print OUTPUT2 "$type,$gc_content\n";
				}
			}
		}
	}
}

%trp_intergenic_hash=();
$trp_intergenic=0;
open REF, "${base_dir}/Data/Reference_files/${species}_reference.fasta";
while(<REF>){
	if(/^([ATGCN]+)/){
		$seq=$1;
		@seq_array=split(//, $seq);
		
		$n=$seq=~tr/N/N/;
		#print "$n\n";
		
		open INPUT, "${base_dir}/Analysis/Core_genome_alignment/${species}_Core_genome_alignment/${species}_core_intergenics.tab";
		while(<INPUT>){
			$line=$_;
			chomp $line;
			if($line =~ /^\S+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)/){
				$beg=$1;$end=$2;$len=$3;$type=$4;
				if($line !~ /^Name\tStart\tEnd\tLength\tType/){
					$indbeg=($beg-1);
					$indend=($end-1);
				
					if($type =~ /Intergenic/){
						$type="Intergenic";
					}
				
					if($len >= 3){
						for($i=$indbeg; $i<=($indend-2); $i++){
							$trp="$seq_array[$i]$seq_array[$i+1]$seq_array[$i+2]";
							$trp_rev=reverse($trp);
						
							@trp_array=($trp, $trp_rev);
							@trp_array=sort(@trp_array);
						
							$trp_intergenic_hash{$trp_array[0]}++;
							$trp_intergenic++;
							#print OUTPUT3 "Intergenic,$trp_array[0]\n";
						}
					}
				
					if($len > 0){
						$gc=0;
						$atgc=0;
						for($i=$indbeg; $i<=$indend; $i++){
							$base=$seq_array[$i];
						
							if($base =~ /G|C/){
								$gc++;
							}
							if($base =~ /A|T|G|C/){
								$atgc++;
							}
						}
					
						$gc_content=($gc / $atgc);
					
						print OUTPUT2 "$type,$gc_content\n";
					}
				}
			}
		}
	}
}

@total_trp_array=keys(%trp_gene_hash);
@total_trp_array=sort(@total_trp_array);

foreach $trp(@total_trp_array){
	$trp_gene_prop=($trp_gene_hash{$trp} / $trp_gene);
	$trp_intergenic_prop=($trp_intergenic_hash{$trp} / $trp_intergenic);
	print OUTPUT "$trp,$trp_gene_prop,$trp_intergenic_prop\n";
}

print "$species base composition summarised.\n";

print LOG "$species base composition summarised.\n";

