#!/usr/bin/bash

species=$1
analysis=$2
base_dir=$3

species_analysis="$species""_""$analysis"

if [ -d "$species_analysis" ]; then
	rm -r "./$species_analysis/"
fi

mkdir "./$species_analysis/"

perl "Variant_summariser.pl" "$species" "$analysis" "$base_dir"

perl "Variant_subsampler.pl" "$species" "$analysis" "$base_dir"

