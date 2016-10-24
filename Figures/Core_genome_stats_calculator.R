library(cowplot)
library(reshape2)

species_array=c("E_coli", "S_aureus", "S_enterica", "S_pneumoniae", "K_pneumoniae", "M_tuberculosis")
species_count=length(species_array)

file_base_a_1="/media/harry/extra/Intergenic_variation_paper/Analysis/Gene_intergenic_coordinates/"
file_base_a_2="_Gene_intergenic_coordinates/"
file_base_a_3_g="_gene_coordinates.tab"
file_base_a_3_i="_intergenic_coordinates.tab"

file_base_b_1="/media/harry/extra/Intergenic_variation_paper/Analysis/Core_genome_alignment/"
file_base_b_2="_Core_genome_alignment/"
file_base_b_3_g="_core_genes.tab"
file_base_b_3_i="_core_intergenics.tab"

species_core_data <- data.frame(Species=character(species_count), Number_genes=numeric(species_count), Number_intergenics=numeric(species_count), Proportion_genes=numeric(species_count), Proportion_intergenic=numeric(species_count), stringsAsFactors=FALSE)

for(i in 1:species_count){
  # if(i == 1){
  species=species_array[i]
  
  g_file=paste(file_base_a_1, species, file_base_a_2, species, file_base_a_3_g, sep="")
  g_data <- read.table(file=g_file, sep="\t", header=TRUE)
  g_data <- g_data[g_data$Type == "CDS", ]
  i_file=paste(file_base_a_1, species, file_base_a_2, species, file_base_a_3_i, sep="")
  i_data <- read.table(file=i_file, sep="\t", header=TRUE)
  
  c_g_file=paste(file_base_b_1, species, file_base_b_2, species, file_base_b_3_g, sep="")
  c_g_data <- read.table(file=c_g_file, sep="\t", header=TRUE)
  c_i_file=paste(file_base_b_1, species, file_base_b_2, species, file_base_b_3_i, sep="")
  c_i_data <- read.table(file=c_i_file, sep="\t", header=TRUE)
  
  g_count <- nrow(g_data)
  i_count <- nrow(i_data)
  
  c_g_count <- nrow(c_g_data)
  c_i_count <- nrow(c_i_data)
  
  c_g_prop <- c_g_count/g_count
  c_i_prop <- c_i_count/i_count
  
  species_core_data$Species[i] <- species
  species_core_data$Number_genes[i] <- g_count
  species_core_data$Number_intergenics[i] <- i_count
  species_core_data$Proportion_genes[i] <- c_g_prop
  species_core_data$Proportion_intergenic[i] <- c_i_prop
  
  #  }
}

species_core_data
