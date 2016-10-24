library(cowplot)
library(reshape2)
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
  
  rows <- nrow(site_data)
  
  Species <- rep(species, rows)
  
  site_data <- cbind(Species, site_data)
  
  species_site_data <- rbind(species_site_data, site_data)
}

species_site_data$Species <- as.character(species_site_data$Species)
species_site_data$Category <- as.character(species_site_data$Category)

species_site_data <- arrange(species_site_data, Species, Category)

Species <- species_site_data$Species
Category <- species_site_data$Category
Count <- species_site_data$Total

species_site_data <- data.frame(Species, Category, Count, stringsAsFactors=FALSE)

species_site_data <- filter(species_site_data, Category == "Intergenic" | Category == "Nonsynonymous" | Category == "Synonymous")

print(species_site_data)

file <- paste(base_dir, "/Figures/Constraint_intergenic.csv", sep="")

constrained_data <- read.csv(file=file, header=TRUE, stringsAsFactors=FALSE)

constrained_data$Synonymous <- 0

constrained_data <- melt(constrained_data, measure.vars=c("Intergenic", "Nonsynonymous", "Synonymous"), id.vars=c("Species"), variable.name="Category", value.name="Constrained_proportion")

constrained_data$Unconstrained_proportion <- 1 - constrained_data$Constrained_proportion

constrained_data$Constrained_sites <- as.integer(constrained_data$Constrained_proportion * species_site_data$Count)
constrained_data$Unconstrained_sites <- as.integer(constrained_data$Unconstrained_proportion * species_site_data$Count)


#constrained_data$Category <- as.character(constrained_data$Category)

#constrained_data$Category[constrained_data$Category == "Intergenic"] <- "Constrained_intergenic"
#constrained_data$Category[constrained_data$Category == "Nonsynonymous"] <- "Constrained_nonsynonymous"

constrained_data <- arrange(constrained_data, Species, Category)

print(constrained_data)

# Summary_data <- read.csv("/media/harry/extra/Intergenic_variation_paper/constraint_for_figure.csv")
# #View(Summary_data)
# 
# Summary_data_long <- melt(Summary_data, id.vars=c("Species", "Order"), measure.vars=c("Synonymous_proportion", "Constrained_N_proportion_total", "Unconstrained_N_proportion_total", "Constrained_I_proportion_total", "Unconstrained_I_proportion_total"), variable.name="Category", value.name="Proportion")
# 
# Summary_data_long_labels <- melt(Summary_data, id.vars="Species", measure.vars=c("Synonymous_sites", "Constrained_N_sites", "Unconstrained_N_sites", "Constrained_I_sites", "Unconstrained_I_sites"), variable.name="Category")
# 
# Summary_data_long <- cbind(Summary_data_long, Summary_data_long_labels$value)
# 
# names(Summary_data_long)[names(Summary_data_long) == 'Summary_data_long_labels$value'] <- 'Sites'
# 
# Summary_data_long$Sites <- as.integer(Summary_data_long$Sites)
# 
# Summary_data_long <- ddply(Summary_data_long, .(Species), transform, pos=cumsum(Proportion) - (0.5 * Proportion))
# 
# print(Summary_data_long)

# species_breaks=c("E_coli", "S_aureus", "S_enterica", "S_pneumoniae", "K_pneumoniae", "M_tuberculosis")
# species_labels=c("E. coli", "S. aureus", "S. enterica", "S. pneumoniae", "K. pneumoniae", "M. tuberculosis")
# 
# category_breaks=c("Synonymous_proportion", "Constrained_N_proportion_total", "Unconstrained_N_proportion_total", "Constrained_I_proportion_total", "Unconstrained_I_proportion_total")
# category_labels=c("Synonymous", "Constrained Non-synonymous", "Unconstrained Non-synonymous", "Constrained Intergenic", "Unconstrained Intergenic")
# 
# Summary_cartoon_plot <- ggplot(Summary_data_long, aes(x=reorder(Species, Order), y=Proportion, fill=Category)) +
#   geom_bar(stat="identity") +
#   geom_text(aes(y=pos, label=Sites)) +
#   scale_x_discrete(breaks=species_breaks, labels=species_labels) +
#   scale_fill_manual(values=c("#339900", "#ff0000", "#ff8080", "#0066ff", "#80b3ff"),
#                     breaks=category_breaks,
#                     labels=category_labels) +
#   labs(x="Species", y="Proportion of genome") +
#   theme(axis.text.x=element_text(face="italic"))
# 
# pdf(file="/media/harry/extra/Intergenic_variation_paper/Figures/Figure_4.pdf", height=10, width=15)
# #tiff(file="/media/harry/extra/Intergenic_variation_paper/Figures/Figure_4.tif", height=10, width=15, units="in", res=100)
# 
# Summary_cartoon_plot
# 
# dev.off()
