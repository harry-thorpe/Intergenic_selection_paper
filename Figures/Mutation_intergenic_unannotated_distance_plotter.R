library(cowplot)
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)
args_count <- length(args)

base_dir <- args[1]

species_array <- NULL

count <- 0
for(i in 2:args_count){
  count <- count + 1
  species_array[count] <- args[i]
}

species_count=length(species_array)

file_base_1="/Analysis/Mutation_intergenic_unannotated_distance/"
file_base_2="_Mutation_intergenic_unannotated_distance/"
file_base_3="_mutation_intergenic_unannotated_distance.tab"

species_mutation_intergenic_unannotated_distance_data <- NULL

for(i in 1:species_count){
  species=species_array[i]
  
  file=paste(base_dir, file_base_1, species, file_base_2, species, file_base_3, sep="")
  
  mutation_intergenic_unannotated_distance_data <- read.csv(file=file, sep="\t", header=TRUE)
  
  rows <- nrow(mutation_intergenic_unannotated_distance_data)
  
  Species <- rep(species, rows)
  
  mutation_intergenic_unannotated_distance_data <- cbind(Species, mutation_intergenic_unannotated_distance_data)
  
  species_mutation_intergenic_unannotated_distance_data <- rbind(species_mutation_intergenic_unannotated_distance_data, mutation_intergenic_unannotated_distance_data)
}

facet_labels=c(E_coli="E. coli", S_aureus="S. aureus", S_enterica="S. enterica", S_pneumoniae="S. pneumoniae", K_pneumoniae="K. pneumoniae", M_tuberculosis="M. tuberculosis")
category_breaks=c("Gene_start_5'", "Gene_end_3'")
category_labels=c("Gene start 5'", "Gene end 3'")

mutation_intergenic_unannotated_distance_plot <- ggplot(species_mutation_intergenic_unannotated_distance_data, aes(x=Distance, y=SNP_density, colour=Category, group=Category)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Species, ncol=3, scales="free", labeller=labeller(Species=facet_labels)) +
  scale_x_continuous(breaks=seq(0, 150, 25)) +
  scale_colour_discrete(breaks=category_breaks, labels=category_labels) +
  labs(x="Distance (bp)", y="SNP density") +
  theme(strip.text.x=element_text(face="italic"))

#out_file_pdf <- paste(base_dir, "/Figures/Figure_6", ".pdf", sep="")
out_file_tif <- paste(base_dir, "/Figures/Figure_6", ".tif", sep="")

#pdf(file=out_file_pdf, height=10, width=15)
tiff(file=out_file_tif, height=10, width=15, units="in", res=100)

mutation_intergenic_unannotated_distance_plot

dev.off()
