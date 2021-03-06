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
library(plyr)

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

file_base_1="/Analysis/Pairwise_dnds_dids_intergenic_annotation/"
file_base_2="_Pairwise_dnds_dids_intergenic_annotation/"
file_base_3="_dnds_dids_intergenic_annotation.csv"

file_base_1_2="/Analysis/Pairwise_dnds_dids/"
file_base_2_2="_Pairwise_dnds_dids/"
file_base_3_2="_dnds_dids.csv"

species_dnds_dids_intergenic_annotation_data <- NULL

for(i in 1:species_count){
  species=species_array[i]
  
  file=paste(base_dir, file_base_1, species, file_base_2, species, file_base_3, sep="")
  
  dnds_dids_intergenic_annotation_data <- read.csv(file=file, header=TRUE)
  
  file=paste(base_dir, file_base_1_2, species, file_base_2_2, species, file_base_3_2, sep="")
  
  dnds_dids_data <- read.csv(file=file, header=TRUE)
  
  dI.dS <- dnds_dids_data$dI.dS
  
  dnds_dids_intergenic_annotation_data <- cbind(dnds_dids_intergenic_annotation_data, dI.dS)
  
  rows <- nrow(dnds_dids_intergenic_annotation_data)
  
  Species <- rep(species, rows)
  
  dnds_dids_intergenic_annotation_data <- cbind(Species, dnds_dids_intergenic_annotation_data)
  
  species_dnds_dids_intergenic_annotation_data <- rbind(species_dnds_dids_intergenic_annotation_data, dnds_dids_intergenic_annotation_data)
}

species_dnds_dids_intergenic_annotation_data_long <- melt(species_dnds_dids_intergenic_annotation_data, measure.vars=c("dN.dS", "dI.dS", "dI.dS_non_coding_RNA", "dI.dS_promoter", "dI.dS_terminator", "dI.dS_unannotated"), variable.name="Category", value.name="dX.dS")

species_dnds_dids_intergenic_annotation_data_long$dS_bin <- cut(species_dnds_dids_intergenic_annotation_data_long$dS, breaks=seq(0, 0.1, 0.0001), labels=seq(0.0001, 0.1, 0.0001))

species_dnds_dids_intergenic_annotation_data_long_summary <- summarySE(species_dnds_dids_intergenic_annotation_data_long, measurevar="dX.dS", groupvars=c("Species", "Category", "dS_bin"))

species_dnds_dids_intergenic_annotation_data_long_summary_summary <- summarySE(species_dnds_dids_intergenic_annotation_data_long_summary, groupvars=c("Species", "Category"), measurevar="dX.dS")

species_dnds_dids_intergenic_annotation_data_long_summary_summary$Constrained_fraction <- (1 - species_dnds_dids_intergenic_annotation_data_long_summary_summary$dX.dS)

####
#species_dnds_dids_intergenic_annotation_data_long_summary_summary$Constrained_fraction <- 1 - species_dnds_dids_intergenic_annotation_data_long_summary_summary$Constrained_fraction
####

species_dnds_dids_intergenic_annotation_data_long_summary_summary_wide <- dcast(species_dnds_dids_intergenic_annotation_data_long_summary_summary, Species ~ Category, value.var="Constrained_fraction")

names(species_dnds_dids_intergenic_annotation_data_long_summary_summary_wide)[names(species_dnds_dids_intergenic_annotation_data_long_summary_summary_wide) == "dN.dS"] <- "Nonsynonymous"
names(species_dnds_dids_intergenic_annotation_data_long_summary_summary_wide)[names(species_dnds_dids_intergenic_annotation_data_long_summary_summary_wide) == "dI.dS"] <- "Intergenic"
names(species_dnds_dids_intergenic_annotation_data_long_summary_summary_wide)[names(species_dnds_dids_intergenic_annotation_data_long_summary_summary_wide) == "dI.dS_non_coding_RNA"] <- "Non_coding_RNA"
names(species_dnds_dids_intergenic_annotation_data_long_summary_summary_wide)[names(species_dnds_dids_intergenic_annotation_data_long_summary_summary_wide) == "dI.dS_promoter"] <- "Promoter"
names(species_dnds_dids_intergenic_annotation_data_long_summary_summary_wide)[names(species_dnds_dids_intergenic_annotation_data_long_summary_summary_wide) == "dI.dS_terminator"] <- "Terminator"
names(species_dnds_dids_intergenic_annotation_data_long_summary_summary_wide)[names(species_dnds_dids_intergenic_annotation_data_long_summary_summary_wide) == "dI.dS_unannotated"] <- "Unannotated"

out_file <- paste(base_dir, "/Figures/Constraint_intergenic.csv", sep="")

write.csv(species_dnds_dids_intergenic_annotation_data_long_summary_summary_wide, file=out_file, quote=FALSE, row.names=FALSE)

dodge <- position_dodge(0.2)

species_breaks=c("E_coli", "S_aureus", "S_enterica", "S_pneumoniae", "K_pneumoniae")
species_labels=c("E. coli", "S. aureus", "S. enterica", "S. pneumoniae", "K. pneumoniae")

category_breaks=c("dN.dS", "dI.dS_non_coding_RNA", "dI.dS_promoter", "dI.dS_terminator", "dI.dS_unannotated", "dI.dS")
category_labels=c("Non-synonymous", "Non coding RNA", "Promoter", "Terminator", "Unannotated", "All intergenic")

Constrained_proportion_intergenic <- ggplot(data=species_dnds_dids_intergenic_annotation_data_long_summary_summary, aes(x=Category, y=Constrained_fraction, fill=Species)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=Constrained_fraction-se, ymax=Constrained_fraction+se), position="dodge") +
  scale_x_discrete(breaks=category_breaks, labels=category_labels) +
  scale_fill_discrete(breaks=species_breaks, labels=species_labels) +
  scale_y_continuous(limits=c(0, 1)) +
  labs(y="Constrained proportion") +
  #labs(y="dX/dS") +
  theme(legend.text=element_text(face="italic"))

#out_file_pdf <- paste(base_dir, "/Figures/Figure_S5", ".pdf", sep="")
out_file_tif <- paste(base_dir, "/Figures/Figure_S5", ".tif", sep="")

#pdf(file=out_file_pdf, height=10, width=15)
tiff(file=out_file_tif, height=10, width=15, units="in", res=100)

Constrained_proportion_intergenic

dev.off()
