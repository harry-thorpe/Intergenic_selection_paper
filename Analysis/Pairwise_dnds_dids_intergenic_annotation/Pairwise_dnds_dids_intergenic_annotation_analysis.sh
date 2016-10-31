#!/usr/bin/bash

species=$1
analysis=$2
base_dir=$3

species_analysis="$species""_""$analysis"

if [ -d "$species_analysis" ]; then
	rm -r "./$species_analysis/"
fi

mkdir "./$species_analysis/"

perl "Alignment_splitter_intergenic_annotation_files.pl" "$species" "$analysis" "$base_dir"

gcc Pairwise_SNP_caller_intergenic_annotation.c -o Pairwise_SNP_caller_intergenic_annotation -lm

category_array=("rbs" "Promoter" "Terminator" "Non_coding_RNA" "Unannotated")
category_file_array=("rbs" "promoter" "terminator" "non_coding_RNA" "unannotated")

category_count=${#category_array[@]}

for ((i=0; i < $category_count ; i++)); do
	
	./Pairwise_SNP_caller_intergenic_annotation "$species" "$analysis" "$base_dir" "${category_array[$i]}" "${category_file_array[$i]}"
done

perl "dnds_dids_combiner_intergenic_annotation.pl" "$species" "$analysis" "$base_dir"

