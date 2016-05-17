#!/usr/bin/bash

##### Start with Intergenic_variation paper, which contains
##### two folders, Analysis, Data, and the Make_paper.sh script.
##### Within the Analysis folder there is a folder for each
##### analysis used in the paper. These folders contain the
##### scripts required for doing each analysis.

base_dir="/media/harry/extra/Intergenic_variation_paper"

cd $base_dir

species_array=("S_aureus" "S_pneumoniae" "E_coli" "S_enterica" "K_pneumoniae" "M_tuberculosis")

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

for species in ${species_array[@]}; do
	for analysis in ${analysis_array[@]}; do
	
		cd "$base_dir/Analysis/$analysis"
	
		bash "$analysis""_analysis.sh" "$species" "$analysis" "$base_dir"
	done
done

