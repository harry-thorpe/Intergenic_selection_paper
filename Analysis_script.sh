#!/usr/bin/bash

# Change this to your own base directory where the code sits.
base_code_dir="/media/harry/extra/Intergenic_selection_paper"

# Change this to your output directory
base_dir="/media/harry/extra/Intergenic_selection_paper_output"

if [ ! -d "$base_dir" ]; then
	mkdir "$base_dir"
fi

#cp -r "$base_code_dir/Data" "$base_dir"
cp -r "$base_code_dir/Analysis" "$base_dir"
cp "$base_code_dir/Analysis_script.sh" "$base_dir"
chmod +x "$base_dir/Analysis/Pairwise_dnds_dids/yn00"
chmod +x "$base_dir/Analysis/Pairwise_dnds_dids_mutation_bias_correction/yn00"
chmod +x "$base_dir/Analysis/Pairwise_dnds_dids_intergenic_annotation_mutation_bias_correction/yn00"

cd "$base_dir"

# Change this to reflect your array of species.
species_array=("S_aureus" "S_pneumoniae" "E_coli" "S_enterica" "K_pneumoniae" "M_tuberculosis")
#species_array=("test")

analysis_array[0]="Gene_intergenic_coordinates"
analysis_array[1]="Core_genome_alignment"
analysis_array[2]="Sequence_summary"
analysis_array[3]="Pairwise_dnds_dids"
analysis_array[4]="Pairwise_dnds_dids_upstream"
analysis_array[5]="Mutation"
analysis_array[6]="Individual_genes_intergenics"
analysis_array[7]="Intergenic_annotation_alignment"
analysis_array[8]="Pairwise_dnds_dids_intergenic_annotation"
analysis_array[9]="Mutation_intergenic_annotation"
analysis_array[10]="Promoter"
analysis_array[11]="Mutation_intergenic_unannotated_distance"
analysis_array[12]="Pairwise_dnds_dids_mutation_bias_correction"
analysis_array[13]="Pairwise_dnds_dids_intergenic_annotation_mutation_bias_correction"

for species in ${species_array[@]}; do
	for analysis in ${analysis_array[@]}; do
	
		cd "$base_dir/Analysis/$analysis"
	
		bash "$analysis""_analysis.sh" "$species" "$analysis" "$base_dir"
	done
done

