library(cowplot)
library(reshape2)
library(plyr)
library(dplyr)

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

file_base_1="/Analysis/Sequence_summary/"
file_base_2="_Sequence_summary/"
file_base_3="_site_counts.csv"

species_site_data <- NULL

for(i in 1:species_count){
  species=species_array[i]
  
  file=paste(base_dir, file_base_1, species, file_base_2, species, file_base_3, sep="")
  
  site_data <- read.csv(file=file, header=TRUE, stringsAsFactors=FALSE)
  
  site_data <- select(site_data, Category, Total)
  
  site_data <- rename(site_data, Count=Total)
  
  site_data <- filter(site_data, Category == "Intergenic" | Category == "Nonsynonymous" | Category == "Synonymous")
  
  Genome_len <- sum(site_data$Count)
  
  site_data$Genome_len <- Genome_len
  
  rows <- nrow(site_data)
  
  Species <- rep(species, rows)
  
  site_data <- cbind(Species, site_data)
  
  species_site_data <- rbind(species_site_data, site_data)
}

species_site_data$Species <- as.character(species_site_data$Species)
species_site_data$Category <- as.character(species_site_data$Category)

species_site_data <- arrange(species_site_data, Species, Category)

#####

file <- paste(base_dir, "/Figures/Constraint_intergenic.csv", sep="")

constrained_data <- read.csv(file=file, header=TRUE, stringsAsFactors=FALSE)

constrained_data$Synonymous <- 0

constrained_data <- melt(constrained_data, measure.vars=c("Intergenic", "Nonsynonymous", "Synonymous"), id.vars=c("Species"), variable.name="Category_1", value.name="Constrained_proportion")

constrained_data <- arrange(constrained_data, Species, Category_1)

constrained_data$Unconstrained_proportion <- 1 - constrained_data$Constrained_proportion

constrained_data$Constrained_sites <- as.integer(constrained_data$Constrained_proportion * species_site_data$Count)
constrained_data$Unconstrained_sites <- as.integer(constrained_data$Unconstrained_proportion * species_site_data$Count)

constrained_data$Constrained_genome_proportion <- (constrained_data$Constrained_sites / species_site_data$Genome_len)
constrained_data$Unconstrained_genome_proportion <- (constrained_data$Unconstrained_sites / species_site_data$Genome_len)

constrained_data_long_1 <- melt(constrained_data, id.vars=c("Species", "Category_1"), measure.vars=c("Constrained_sites", "Unconstrained_sites"), variable.name="Category_2", value.name="Sites")

constrained_data_long_1$Category_2 <- as.character(constrained_data_long_1$Category_2)

constrained_data_long_1$Category_2[constrained_data_long_1$Category_2 == "Constrained_sites"] <- "Constrained"
constrained_data_long_1$Category_2[constrained_data_long_1$Category_2 == "Unconstrained_sites"] <- "Unconstrained"

constrained_data_long_2 <- melt(constrained_data, id.vars=c("Species", "Category_1"), measure.vars=c("Constrained_genome_proportion", "Unconstrained_genome_proportion"), variable.name="Category_2", value.name="Genome_proportion")

constrained_data_long_2$Category_2 <- as.character(constrained_data_long_2$Category_2)

constrained_data_long_2$Category_2[constrained_data_long_2$Category_2 == "Constrained_genome_proportion"] <- "Constrained"
constrained_data_long_2$Category_2[constrained_data_long_2$Category_2 == "Unconstrained_genome_proportion"] <- "Unconstrained"

constrained_data_long <- merge(constrained_data_long_1, constrained_data_long_2)

constrained_data_long <- constrained_data_long[! constrained_data_long$Sites == 0, ]

constrained_data_long$Category <- paste(constrained_data_long$Category_1, constrained_data_long$Category_2, sep="_")

constrained_data_long$Order <- NULL

constrained_data_long$Order[constrained_data_long$Category == "Synonymous_Unconstrained"] <- 1
constrained_data_long$Order[constrained_data_long$Category == "Nonsynonymous_Constrained"] <- 2
constrained_data_long$Order[constrained_data_long$Category == "Nonsynonymous_Unconstrained"] <- 3
constrained_data_long$Order[constrained_data_long$Category == "Intergenic_Constrained"] <- 4
constrained_data_long$Order[constrained_data_long$Category == "Intergenic_Unconstrained"] <- 5

constrained_data_long <- arrange(constrained_data_long, Species, Order)

print(constrained_data_long)

constrained_data_long <- ddply(constrained_data_long, .(Species), transform, pos=cumsum(Genome_proportion) - (0.5 * Genome_proportion))

print(constrained_data_long)

species_breaks=c("E_coli", "S_aureus", "S_enterica", "S_pneumoniae", "K_pneumoniae")
species_labels=c("E. coli", "S. aureus", "S. enterica", "S. pneumoniae", "K. pneumoniae")

category_values=c("Synonymous_Unconstrained"="#339900", "Nonsynonymous_Constrained"="#ff0000", "Nonsynonymous_Unconstrained"="#ff8080", "Intergenic_Constrained"="#0066ff", "Intergenic_Unconstrained"="#80b3ff")
#category_breaks=c("Synonymous_Unconstrained", "Nonsynonymous_Constrained", "Nonsynonymous_Unconstrained", "Intergenic_Constrained", "Intergenic_Unconstrained")
category_labels=c("Synonymous_Unconstrained"="Synonymous", "Nonsynonymous_Constrained"="Non-synonymous constrained", "Nonsynonymous_Unconstrained"="Non-synonymous unconstrained", "Intergenic_Constrained"="Intergenic constrained", "Intergenic_Unconstrained"="Intergenic unconstrained")

Summary_cartoon_plot <- ggplot(constrained_data_long, aes(x=Species, y=Genome_proportion, fill=Category)) +
  geom_bar(stat="identity") +
  geom_text(aes(y=pos, label=Sites)) +
  scale_x_discrete(breaks=species_breaks, labels=species_labels) +
  scale_fill_manual(values=category_values,
                    #breaks=category_breaks,
                    labels=category_labels) +
  labs(x="Species", y="Proportion of genome") +
  theme(axis.text.x=element_text(face="italic"))

#out_file_pdf <- paste(base_dir, "/Figures/Figure_4", ".pdf", sep="")
out_file_tif <- paste(base_dir, "/Figures/Figure_4", ".tif", sep="")

#pdf(file=out_file_pdf, height=10, width=15)
tiff(file=out_file_tif, height=10, width=15, units="in", res=100)

Summary_cartoon_plot

dev.off()
