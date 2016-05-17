#!/usr/bin/bash

species=$1
analysis=$2
base_dir=$3

species_analysis="$species""_""$analysis"

rm -r "./$species_analysis/"

mkdir "./$species_analysis/"

perl "Alignment_creator_core_intergenic_annotations.pl" "$species" "$analysis" "$base_dir"

