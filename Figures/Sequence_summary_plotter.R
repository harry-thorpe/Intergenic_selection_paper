## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#####
library(cowplot)

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
file_base_3="_GC_content_summary.csv"

species_GC_data <- NULL

for(i in 1:species_count){
  species=species_array[i]
  
  file=paste(base_dir, file_base_1, species, file_base_2, species, file_base_3, sep="")
  
  GC_data <- read.csv(file=file, header=TRUE)
  
  rows <- nrow(GC_data)
  
  Species <- rep(species, rows)
  
  GC_data <- cbind(Species, GC_data)
  
  species_GC_data <- rbind(species_GC_data, GC_data)
}

species_GC_data_summary <- summarySE(species_GC_data, groupvars=c("Species", "Category"), measurevar="GC_content")

species_breaks=c("E_coli", "S_aureus", "S_enterica", "S_pneumoniae", "K_pneumoniae", "M_tuberculosis")
species_labels=c("E. coli", "S. aureus", "S. enterica", "S. pneumoniae", "K. pneumoniae", "M. tuberculosis")

category_breaks=c("Gene", "Intergenic")
category_labels=c("Gene\n", "Intergenic\n")

gc_plot <- ggplot() +
  geom_point(data=species_GC_data_summary, aes(x=Category, y=GC_content, colour=Species)) +
  geom_errorbar(data=species_GC_data_summary, aes(x=Category, ymin=GC_content-se, ymax=GC_content+se, width=0.1, colour=Species)) +
  geom_line(data=species_GC_data_summary, aes(x=Category, y=GC_content, colour=Species, group=Species)) +
  scale_y_continuous(limits=c(0, 0.7), breaks=seq(0, 0.7, 0.1)) +
  scale_x_discrete(breaks=category_breaks, labels=category_labels) +
  scale_colour_discrete(breaks=species_breaks, labels=species_labels) +
  labs(y="GC content") +
  theme(legend.text=element_text(face="italic"))

#####

file_base_1="/Analysis/Core_genome_alignment/"
file_base_2="_Core_genome_alignment/"
file_base_3="_core_intergenics.tab"

species_intergenic_data <- NULL

for(i in 1:species_count){
  species=species_array[i]
  
  file=paste(base_dir, file_base_1, species, file_base_2, species, file_base_3, sep="")
  
  intergenic_data <- read.csv(file=file, sep="\t", header=TRUE)
  
  rows <- nrow(intergenic_data)
  
  Species <- rep(species, rows)
  
  intergenic_data <- cbind(Species, intergenic_data)
  
  species_intergenic_data <- rbind(species_intergenic_data, intergenic_data)
}

species_intergenic_data$Category <- species_intergenic_data$Type

species_intergenic_data$Category <- as.character(species_intergenic_data$Category)

for(i in 1:nrow(species_intergenic_data)){
  if(species_intergenic_data$Type[i] == "Intergenic_co-oriented_F"){
    species_intergenic_data$Category[i] <- "Intergenic_CO"
  }else if(species_intergenic_data$Type[i] == "Intergenic_co-oriented_R"){
    species_intergenic_data$Category[i] <- "Intergenic_CO"
  }else if(species_intergenic_data$Type[i] == "Intergenic_double_promoter"){
    species_intergenic_data$Category[i] <- "Intergenic_DP"
  }else if(species_intergenic_data$Type[i] == "Intergenic_double_terminator"){
    species_intergenic_data$Category[i] <- "Intergenic_DT"
  }
}

species_intergenic_data_summary <- summarySE(species_intergenic_data, groupvars=c("Species", "Category"), measurevar="Length")

#####

file_base_1="/Analysis/Core_genome_alignment/"
file_base_2="_Core_genome_alignment/"
file_base_3="_core_genes.tab"

species_gene_data <- NULL

for(i in 1:species_count){
  species=species_array[i]
  
  file=paste(base_dir, file_base_1, species, file_base_2, species, file_base_3, sep="")
  
  gene_data <- read.csv(file=file, sep="\t", header=TRUE)
  
  rows <- nrow(gene_data)
  
  Species <- rep(species, rows)
  
  gene_data <- cbind(Species, gene_data)
  
  species_gene_data <- rbind(species_gene_data, gene_data)
}

species_gene_data$Category <- species_gene_data$Type

species_gene_data$Category <- as.character(species_gene_data$Category)

for(i in 1:nrow(species_gene_data)){
  if(species_gene_data$Type[i] == "CDS"){
    species_gene_data$Category[i] <- "Gene"
  }
}

species_gene_data_summary <- summarySE(species_gene_data, groupvars=c("Species", "Category"), measurevar="Length")

#####

species_length_data_summary <- rbind(species_gene_data_summary, species_intergenic_data_summary)

category_breaks=c("Gene", "Intergenic_CO", "Intergenic_DP", "Intergenic_DT")
category_labels=c("Gene", "Intergenic\nCO", "Intergenic\nDP", "Intergenic\nDT")

dodge <- position_dodge(0.6)

length_plot <- ggplot() +
  geom_boxplot(data=species_length_data_summary, aes(x=Category, y=Length), outlier.size=NA) +
  geom_point(data=species_length_data_summary, aes(x=Category, y=Length, colour=Species), position=dodge) +
  geom_errorbar(data=species_length_data_summary, aes(x=Category, ymin=Length-se, ymax=Length+se, width=0.4, colour=Species), position=dodge) +
  scale_y_continuous(limits=c(0, 1000), breaks=seq(0, 1000, 200)) +
  scale_x_discrete(breaks=category_breaks, labels=category_labels) +
  labs(y="Length (bp)") +
  theme(legend.position="none")

#####

file_base_1="/Analysis/Sequence_summary/"
file_base_2="_Sequence_summary/"
file_base_3="_site_counts.csv"

species_site_data <- NULL

for(i in 1:species_count){
  species=species_array[i]
  
  file=paste(base_dir, file_base_1, species, file_base_2, species, file_base_3, sep="")
  
  site_data <- read.csv(file=file, header=TRUE)
  
  site_data <- site_data[site_data$Category != "Nonsense", ]
  names(site_data)[names(site_data) == "Category"] <- "Site"
  
  site_data$Total_GC_content <- ((sum(site_data$G) + sum(site_data$C)) / sum(site_data$Total))
  
  rows <- nrow(site_data)
  
  Species <- rep(species, rows)
  
  site_data <- cbind(Species, site_data)
  
  species_site_data <- rbind(species_site_data, site_data)
}

site_labels=c("Intergenic", "Non-synonymous", "Synonymous")
site_breaks=c("Intergenic", "Nonsynonymous", "Synonymous")

species_site_plot <- ggplot(species_site_data, aes(x=Total_GC_content, y=GC_content, colour=Site)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  scale_x_continuous(limits=c(0.3,0.7), breaks=seq(0.3, 0.7, 0.1)) +
  scale_y_continuous(limits=c(0,0.9), breaks=seq(0, 0.9, 0.1)) +
  scale_colour_discrete(breaks=site_breaks, labels=site_labels) +
  #coord_equal(ratio=0.7) +
  labs(x="Genome GC content", y="Site GC content")

#####

sequence_summary_plot <- ggdraw() +
  draw_plot(length_plot, 0, 0, 0.3, 1) +
  draw_plot(gc_plot, 0.3, 0, 0.3, 1) +
  draw_plot(species_site_plot, 0.6, 0, 0.4, 1) +
  draw_plot_label(c("a", "b", "c"), c(0, 0.3, 0.6), c(1, 1, 1))

#out_file_pdf <- paste(base_dir, "/Figures/Figure_1", ".pdf", sep="")
out_file_tif <- paste(base_dir, "/Figures/Figure_1", ".tif", sep="")

#pdf(file=out_file_pdf, height=7, width=15)
tiff(file=out_file_tif, height=7, width=15, units="in", res=100)

sequence_summary_plot

dev.off()
