#!/usr/bin/bash

species=$1
analysis=$2
base_dir=$3

species_analysis="$species""_""$analysis"

if [ -d "$species_analysis" ]; then
	rm -r "./$species_analysis/"
fi

mkdir "./$species_analysis/"

perl "Simulator.pl" "$species" "$analysis" "$base_dir"

# Calculate dN/dS
mkdir "$base_dir/Analysis/$analysis/$species_analysis/dnds_tmp"

cd "$base_dir/Analysis/$analysis/$species_analysis/dnds_tmp"

cp "$base_dir/Analysis/$analysis/yn00.ctl" "$base_dir/Analysis/$analysis/$species_analysis/dnds_tmp"

cp "$base_dir/Analysis/$analysis/$species_analysis/${species}_core_gene_alignment_simulated.fasta" "$base_dir/Analysis/$analysis/$species_analysis/dnds_tmp"

mv "$base_dir/Analysis/$analysis/$species_analysis/dnds_tmp/${species}_core_gene_alignment_simulated.fasta" "$base_dir/Analysis/$analysis/$species_analysis/dnds_tmp/ali.fasta"

${base_dir}/Analysis/$analysis/yn00 yn00.ctl

cd "$base_dir/Analysis/$analysis"

perl "dnds_parser.pl" "$species" "$analysis" "$base_dir"

rm -r "$base_dir/Analysis/$analysis/$species_analysis/dnds_tmp"

# Calculate dI
perl "Alignment_splitter_intergenic_upstream_files.pl" "$species" "$analysis" "$base_dir"

gcc Pairwise_SNP_caller_intergenic_upstream.c -o Pairwise_SNP_caller_intergenic_upstream -lm

category_array=("all" "upstream_30")
category_file_array=("all" "upstream_30")

category_count=${#category_array[@]}

for ((i=0; i < $category_count ; i++)); do
	
	./Pairwise_SNP_caller_intergenic_upstream "$species" "$analysis" "$base_dir" "${category_array[$i]}" "${category_file_array[$i]}"
done

# Combine dN/dS and dI
perl "dnds_dids_combiner.pl" "$species" "$analysis" "$base_dir"

