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

file_base_1="/Analysis/Individual_genes_intergenics/"
file_base_2="_Individual_genes_intergenics/"
file_base_3="_core_gene_intergenic_PSM.tab"

species_core_gene_intergenic_PSM <- NULL

for(i in 1:species_count){
  species=species_array[i]
  
  file=paste(base_dir, file_base_1, species, file_base_2, species, file_base_3, sep="")
  
  core_gene_intergenic_PSM <- read.csv(file=file, sep="\t", header=TRUE)
  
  rows <- nrow(core_gene_intergenic_PSM)
  
  Species <- rep(species, rows)
  
  core_gene_intergenic_PSM <- cbind(Species, core_gene_intergenic_PSM)
  
  species_core_gene_intergenic_PSM <- rbind(species_core_gene_intergenic_PSM, core_gene_intergenic_PSM)
}

species_core_gene_intergenic_PSM$Order <- species_core_gene_intergenic_PSM$PSM

row_count <- nrow(species_core_gene_intergenic_PSM)

for(i in 1:row_count){
  if(species_core_gene_intergenic_PSM$Category[i] == "Synonymous"){
    species_core_gene_intergenic_PSM$Order[i] <- 1
  }else if(species_core_gene_intergenic_PSM$Category[i] == "Intergenic"){
    species_core_gene_intergenic_PSM$Order[i] <- 2
  }else if(species_core_gene_intergenic_PSM$Category[i] == "Nonsynonymous"){
    species_core_gene_intergenic_PSM$Order[i] <- 3
  }else if(species_core_gene_intergenic_PSM$Category[i] == "Nonsense"){
    species_core_gene_intergenic_PSM$Order[i] <- 4
  }
}

facet_labels=c(E_coli="E. coli", S_aureus="S. aureus", S_enterica="S. enterica", S_pneumoniae="S. pneumoniae", K_pneumoniae="K. pneumoniae", M_tuberculosis="M. tuberculosis")
category_breaks=c("Synonymous", "Intergenic", "Nonsynonymous", "Nonsense")
category_labels=c("Synonymous", "Intergenic", "Non-synonymous", "Nonsense")

core_gene_intergenic_PSM_plot <- ggplot(species_core_gene_intergenic_PSM, aes(x=reorder(Category, Order), y=PSM, colour=Category)) +
  geom_boxplot(notch=TRUE) +
  facet_wrap(~Species, ncol=3, labeller=labeller(Species=facet_labels)) +
  scale_x_discrete(breaks=category_breaks, labels=category_labels) +
  labs(x="Category", y="PSM") +
  theme(strip.text.x=element_text(face="italic"),
        legend.position="none")

#out_file_pdf <- paste(base_dir, "/Figures/Figure_S1", ".pdf", sep="")
out_file_tif <- paste(base_dir, "/Figures/Figure_S1", ".tif", sep="")

#pdf(file=out_file_pdf, height=10, width=15)
tiff(file=out_file_tif, height=10, width=15, units="in", res=100)

core_gene_intergenic_PSM_plot

dev.off()
