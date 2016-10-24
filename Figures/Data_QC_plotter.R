library(cowplot)
library(reshape2)

species_array=c("E_coli", "S_aureus", "S_enterica", "S_pneumoniae", "K_pneumoniae", "M_tuberculosis")
species_count=length(species_array)

file_base_1="/media/harry/extra/Intergenic_variation_paper/Analysis/Data_QC/"
file_base_2="_Data_QC/"
file_base_3="_variants_subsampled.tab"

species_data_qc <- NULL

for(i in 1:species_count){
  # if(i == 1){
  species=species_array[i]
  
  file=paste(file_base_1, species, file_base_2, species, file_base_3, sep="")
  
  data_qc <- read.delim(file, header=TRUE)
  
  rows <- nrow(data_qc)
  
  Species <- rep(species, rows)
  
  data_qc <- cbind(Species, data_qc)
  
  species_data_qc <- rbind(species_data_qc, data_qc)
  # }
}

species_data_qc$Category <- "Non_singleton"
species_data_qc$Category[species_data_qc$Count == 1] <- "Singleton"
#species_data_qc$Category[species_data_qc$Count == 2] <- "Doubleton"
#species_data_qc$Category[species_data_qc$Count == 3] <- "Tripleton"

ggplot(species_data_qc, aes(x=Depth, colour=Category)) +
  geom_density() +
  scale_x_continuous(limits=c(0,200)) +
  facet_wrap(~Species, ncol=3, scales="free")

ggplot(species_data_qc, aes(x=Depth, fill=Category)) +
  geom_histogram() +
  scale_x_continuous(limits=c(0,200)) +
  facet_wrap(~Species, ncol=3, scales="free")


ggplot(species_data_qc, aes(x=Support, colour=Category)) +
  geom_density() +
  scale_x_continuous(limits=c(0.9,1.01)) +
  facet_wrap(~Species, ncol=3, scales="free")

ggplot(species_data_qc, aes(x=Support, fill=Category)) +
  geom_histogram() +
  facet_wrap(~Species, ncol=3, scales="free")


ggplot(species_data_qc, aes(x=Score, colour=Category)) +
  geom_density() +
  scale_x_continuous(limits=c(200,225)) +
  facet_wrap(~Species, ncol=3, scales="free")

ggplot(species_data_qc, aes(x=Score, fill=Category)) +
  geom_histogram() +
  facet_wrap(~Species, ncol=3, scales="free")
