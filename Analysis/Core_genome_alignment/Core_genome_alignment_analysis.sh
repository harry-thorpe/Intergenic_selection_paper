#!/usr/bin/bash

species=$1
analysis=$2
base_dir=$3

species_analysis="$species""_""$analysis"

rm -r "./$species_analysis/"

mkdir "./$species_analysis/"

perl "Alignment_splitter_genes.pl" "$species" "$analysis" "$base_dir"

perl "Alignment_splitter_intergenics.pl" "$species" "$analysis" "$base_dir"

perl "Alignment_checker_core_genes.pl" "$species" "$analysis" "$base_dir"

perl "Alignment_checker_core_intergenics.pl" "$species" "$analysis" "$base_dir"

rm -r "./$species_analysis/${species}""_gene_files"

rm -r "./$species_analysis/${species}""_intergenic_files"

perl "Alignment_creator_core_genes.pl" "$species" "$analysis" "$base_dir"

perl "Alignment_creator_core_intergenics.pl" "$species" "$analysis" "$base_dir"

