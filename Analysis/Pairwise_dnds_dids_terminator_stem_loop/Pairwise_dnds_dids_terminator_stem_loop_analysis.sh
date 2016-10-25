#!/usr/bin/bash

species=$1
analysis=$2
base_dir=$3

species_analysis="$species""_""$analysis"

if [ -d "$species_analysis" ]; then
	rm -r "./$species_analysis/"
fi

mkdir "./$species_analysis/"

perl "Alignment_creator_terminator_stem_loop.pl" "$species" "$analysis" "$base_dir"

perl "Alignment_splitter_terminator_stem_loop_files.pl" "$species" "$analysis" "$base_dir"

gcc Pairwise_SNP_caller_terminator_stem_loop.c -o Pairwise_SNP_caller_terminator_stem_loop -lm

category_array=("terminator_stem" "terminator_loop")
category_file_array=("terminator_stem" "terminator_loop")

category_count=${#category_array[@]}

for ((i=0; i < $category_count ; i++)); do
	
	./Pairwise_SNP_caller_terminator_stem_loop "$species" "$analysis" "$base_dir" "${category_array[$i]}" "${category_file_array[$i]}"
done

perl "dnds_dids_combiner_terminator_stem_loop.pl" "$species" "$analysis" "$base_dir"

