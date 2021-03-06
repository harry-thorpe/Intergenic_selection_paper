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

${base_dir}/Analysis/$analysis/yn00 yn00.ctl &> /dev/null

cd "$base_dir/Analysis/$analysis"

perl "dnds_parser.pl" "$species" "$analysis" "$base_dir"

#rm -r "$base_dir/Analysis/$analysis/$species_analysis/dnds_tmp"

# Calculate dI
perl "Alignment_splitter_intergenic_files.pl" "$species" "$analysis" "$base_dir"

gcc Pairwise_SNP_caller_intergenic.c -o Pairwise_SNP_caller_intergenic -lm

./Pairwise_SNP_caller_intergenic "$species" "$analysis" "$base_dir"

# Combine dN/dS and dI
perl "dnds_dids_combiner.pl" "$species" "$analysis" "$base_dir"

