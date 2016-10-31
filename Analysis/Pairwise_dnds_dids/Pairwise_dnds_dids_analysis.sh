#!/usr/bin/bash

species=$1
analysis=$2
base_dir=$3

species_analysis="$species""_""$analysis"

if [ -d "$species_analysis" ]; then
	rm -r "./$species_analysis/"
fi

mkdir "./$species_analysis/"

threshold_array=("0" "95" "99")

for threshold in ${threshold_array[@]}; do
	
	if [ "$threshold" -ne "95" ]; then
		mkdir "${base_dir}/Analysis/${analysis}/${species}_${analysis}/threshold_$threshold"
	fi

	perl "STO_remover_genes.pl" "$species" "$analysis" "$base_dir" "$threshold"

	perl "Shuffler_genes.pl" "$species" "$analysis" "$base_dir" "$threshold"

	perl "Alignment_splitter_gene_files.pl" "$species" "$analysis" "$base_dir" "$threshold"

	perl "Alignment_splitter_intergenic_files.pl" "$species" "$analysis" "$base_dir" "$threshold"

	parallel --no-notice "bash yn00_executor.sh $species $analysis $base_dir $threshold {}" ::: {a..b}

	gcc Pairwise_SNP_caller_intergenic.c -o Pairwise_SNP_caller_intergenic -lm

	./Pairwise_SNP_caller_intergenic "$species" "$analysis" "$base_dir" "$threshold"

	perl "dnds_dids_combiner.pl" "$species" "$analysis" "$base_dir" "$threshold"
done

