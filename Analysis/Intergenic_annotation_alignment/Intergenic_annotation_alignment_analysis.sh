#!/usr/bin/bash

species=$1
analysis=$2
base_dir=$3

species_analysis="$species""_""$analysis"

if [ -d "$species_analysis" ]; then
	rm -r "./$species_analysis/"
fi

mkdir "./$species_analysis/"

perl "rbs_finder_input_creator.pl" "$species" "$analysis" "$base_dir"

perl "rbs_finder.pl" "$base_dir/Data/Reference_files/${species}_reference.fasta" "$base_dir/Analysis/$analysis/$species_analysis/${species}_gene_coordinates_for_rbs_finder.tab" "$base_dir/Analysis/$analysis/$species_analysis/${species}_rbs_finder_output.tab" "15"

perl "rbs_finder_output_parser.pl" "$species" "$analysis" "$base_dir"

perl "Alignment_creator_core_intergenic_annotations.pl" "$species" "$analysis" "$base_dir"

