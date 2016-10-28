#!/usr/bin/bash

species=$1
analysis=$2
base_dir=$3

species_analysis="$species""_""$analysis"

if [ -d "$species_analysis" ]; then
	rm -r "./$species_analysis/"
fi

mkdir "./$species_analysis/"

perl "Alignment_splitter_genes.pl" "$species" "$analysis" "$base_dir"

perl "Alignment_splitter_intergenics.pl" "$species" "$analysis" "$base_dir"

threshold_array=("0" "95" "99")

for threshold in ${threshold_array[@]}; do
	
	if [ "$threshold" -ne "95" ]; then
		mkdir "${base_dir}/Analysis/${analysis}/${species}_${analysis}/threshold_$threshold"
	fi

	perl "Alignment_checker_core_genes.pl" "$species" "$analysis" "$base_dir" "$threshold"

	perl "Alignment_checker_core_intergenics.pl" "$species" "$analysis" "$base_dir" "$threshold"

	perl "Alignment_creator_core_genes.pl" "$species" "$analysis" "$base_dir" "$threshold"

	perl "Alignment_creator_core_intergenics.pl" "$species" "$analysis" "$base_dir" "$threshold"
done

