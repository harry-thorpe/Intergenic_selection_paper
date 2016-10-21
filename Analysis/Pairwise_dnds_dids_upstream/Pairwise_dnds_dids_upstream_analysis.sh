#!/usr/bin/bash

species=$1
analysis=$2
base_dir=$3

species_analysis="$species""_""$analysis"

if [ -d "$species_analysis" ]; then
	rm -r "./$species_analysis/"
fi

mkdir "./$species_analysis/"

upstream_len_array[0]=30

for i in ${upstream_len_array[@]}; do
	perl "Alignment_creator_core_intergenics_upstream.pl" "$species" "$analysis" "$base_dir" "$i"

	perl "Alignment_splitter_intergenic_upstream_files.pl" "$species" "$analysis" "$base_dir" "$i"

	gcc Pairwise_SNP_caller_intergenic_upstream.c -o Pairwise_SNP_caller_intergenic_upstream -lm

	./Pairwise_SNP_caller_intergenic_upstream "$species" "$analysis" "$base_dir" "$i"
done

perl "dnds_dids_combiner_upstream.pl" "$species" "$analysis" "$base_dir"

