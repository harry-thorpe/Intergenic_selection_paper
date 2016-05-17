#!/usr/bin/bash

species=$1
analysis=$2
base_dir=$3

species_analysis="$species""_""$analysis"

rm -r "./$species_analysis/"

mkdir "./$species_analysis/"

perl "Individual_genes_intergenics_PSM_calculator.pl" "$species" "$analysis" "$base_dir"

