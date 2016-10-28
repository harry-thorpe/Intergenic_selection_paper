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

file_base_1="/Analysis/Mutation/"
file_base_2="_Mutation/"
file_base_3="_mutation_bias.tab"

species_mutation_bias_data <- NULL

for(i in 1:species_count){
  species=species_array[i]
  
  file=paste(base_dir, file_base_1, species, file_base_2, species, file_base_3, sep="")
  
  mutation_bias_data <- read.delim(file=file, header=TRUE)
  
  rows <- nrow(mutation_bias_data)
  
  Species <- rep(species, rows)
  
  mutation_bias_data <- cbind(Species, mutation_bias_data)
  
  species_mutation_bias_data <- rbind(species_mutation_bias_data, mutation_bias_data)
}

#species_mutation_bias_data <- species_mutation_bias_data[species_mutation_bias_data$Category != "Nonsense", ]
species_mutation_bias_data <- species_mutation_bias_data[species_mutation_bias_data$Category == "All", ]

mutation_breaks=c("AC", "AG", "AT", "CA", "CG", "CT")
#mutation_labels=c("A/T -> C/G", "A/T -> G/C", "A/T -> T/A", "C/G -> A/T", "C/G -> G/C", "C/G -> T/A")
mutation_labels=c("A->C\nor\nT->G", "A->G\nor\nT->C", "A->T\nor\nT->A", "C->A\nor\nG->T", "C->G\nor\nG->C", "C->T\nor\nG->A")

facet_labels=c(E_coli="E. coli", S_aureus="S. aureus", S_enterica="S. enterica", S_pneumoniae="S. pneumoniae", K_pneumoniae="K. pneumoniae", M_tuberculosis="M. tuberculosis")

mutation_bias_plot <- ggplot(species_mutation_bias_data, aes(x=Mutation, y=Count, fill=Category)) +
  geom_bar(stat="identity", position="dodge") +
  scale_x_discrete(breaks=mutation_breaks, labels=mutation_labels) +
  facet_wrap(~Species, ncol=3, labeller=labeller(Species=facet_labels)) +
  theme(strip.text.x=element_text(face="italic"))

#out_file_pdf <- paste(base_dir, "/Figures/Figure_S9", ".pdf", sep="")
out_file_tif <- paste(base_dir, "/Figures/Figure_S9", ".tif", sep="")

#pdf(file=out_file_pdf, height=10, width=15)
tiff(file=out_file_tif, height=10, width=15, units="in", res=100)

mutation_bias_plot

dev.off()
