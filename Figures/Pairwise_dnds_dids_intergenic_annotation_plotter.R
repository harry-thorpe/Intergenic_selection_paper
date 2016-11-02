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

file_base_1="/Analysis/Pairwise_dnds_dids_intergenic_annotation/"
file_base_2="_Pairwise_dnds_dids_intergenic_annotation/"
file_base_3="_dnds_dids_intergenic_annotation.csv"

file_base_1_1="/Analysis/Pairwise_dnds_dids_intergenic_annotation_mutation_bias_correction/"
file_base_2_2="_Pairwise_dnds_dids_intergenic_annotation_mutation_bias_correction/"
file_base_3_3="_dnds_dids_intergenic_annotation_simulated.csv"

species_dnds_dids_intergenic_annotation_data <- NULL

for(i in 1:species_count){
  species=species_array[i]
  
  file=paste(base_dir, file_base_1_1, species, file_base_2_2, species, file_base_3_3, sep="")
  
  dnds_dids_intergenic_annotation_data <- read.csv(file=file, header=TRUE)
  
  mean_dN.dS <- mean(dnds_dids_intergenic_annotation_data$dN.dS)
  mean_dI.dS_rbs <- mean(dnds_dids_intergenic_annotation_data$dI.dS_rbs)
  mean_dI.dS_non_coding_RNA <- mean(dnds_dids_intergenic_annotation_data$dI.dS_non_coding_RNA)
  mean_dI.dS_promoter <- mean(dnds_dids_intergenic_annotation_data$dI.dS_promoter)
  mean_dI.dS_terminator <- mean(dnds_dids_intergenic_annotation_data$dI.dS_terminator)
  mean_dI.dS_unannotated <- mean(dnds_dids_intergenic_annotation_data$dI.dS_unannotated)
  
  file=paste(base_dir, file_base_1, species, file_base_2, species, file_base_3, sep="")
  
  dnds_dids_intergenic_annotation_data <- read.csv(file=file, header=TRUE)
  
  dnds_dids_intergenic_annotation_data$dN.dS <- dnds_dids_intergenic_annotation_data$dN.dS/mean_dN.dS
  dnds_dids_intergenic_annotation_data$dI.dS_rbs <- dnds_dids_intergenic_annotation_data$dI.dS_rbs/mean_dI.dS_rbs
  dnds_dids_intergenic_annotation_data$dI.dS_non_coding_RNA <- dnds_dids_intergenic_annotation_data$dI.dS_non_coding_RNA/mean_dI.dS_non_coding_RNA
  dnds_dids_intergenic_annotation_data$dI.dS_promoter <- dnds_dids_intergenic_annotation_data$dI.dS_promoter/mean_dI.dS_promoter
  dnds_dids_intergenic_annotation_data$dI.dS_terminator <- dnds_dids_intergenic_annotation_data$dI.dS_terminator/mean_dI.dS_terminator
  dnds_dids_intergenic_annotation_data$dI.dS_unannotated <- dnds_dids_intergenic_annotation_data$dI.dS_unannotated/mean_dI.dS_unannotated
  
  rows <- nrow(dnds_dids_intergenic_annotation_data)
  
  Species <- rep(species, rows)
  
  dnds_dids_intergenic_annotation_data <- cbind(Species, dnds_dids_intergenic_annotation_data)
  
  species_dnds_dids_intergenic_annotation_data <- rbind(species_dnds_dids_intergenic_annotation_data, dnds_dids_intergenic_annotation_data)
}

species_dnds_dids_intergenic_annotation_data_long <- melt(species_dnds_dids_intergenic_annotation_data, measure.vars=c("dN.dS", "dI.dS_rbs", "dI.dS_non_coding_RNA", "dI.dS_promoter", "dI.dS_terminator", "dI.dS_unannotated"), variable.name="Category", value.name="dX.dS")

species_dnds_dids_intergenic_annotation_data_long$dS_bin <- cut(species_dnds_dids_intergenic_annotation_data_long$dS, breaks=seq(0, 0.1, 0.0001), labels=seq(0.0001, 0.1, 0.0001))

species_dnds_dids_intergenic_annotation_data_long$dS_bin <- as.numeric(as.character(species_dnds_dids_intergenic_annotation_data_long$dS_bin))

species_dnds_dids_intergenic_annotation_data_long_summary <- summarySE(species_dnds_dids_intergenic_annotation_data_long, measurevar="dX.dS", groupvars=c("Species", "Category", "dS_bin"))

species_dnds_dids_intergenic_annotation_data_long_summary_wide <- dcast(species_dnds_dids_intergenic_annotation_data_long_summary, Category + dS_bin ~ Species, value.var="dX.dS")

species_dnds_dids_intergenic_annotation_data_long_summary_wide_E_coli <- dcast(species_dnds_dids_intergenic_annotation_data_long_summary_wide, dS_bin ~ Category, value.var="E_coli")

species_dnds_dids_intergenic_annotation_data_long_summary_wide_S_aureus <- dcast(species_dnds_dids_intergenic_annotation_data_long_summary_wide, dS_bin ~ Category, value.var="S_aureus")

species_dnds_dids_intergenic_annotation_data_long_summary_wide_S_enterica <- dcast(species_dnds_dids_intergenic_annotation_data_long_summary_wide, dS_bin ~ Category, value.var="S_enterica")

species_dnds_dids_intergenic_annotation_data_long_summary_wide_S_pneumoniae <- dcast(species_dnds_dids_intergenic_annotation_data_long_summary_wide, dS_bin ~ Category, value.var="S_pneumoniae")

species_dnds_dids_intergenic_annotation_data_long_summary_wide_K_pneumoniae <- dcast(species_dnds_dids_intergenic_annotation_data_long_summary_wide, dS_bin ~ Category, value.var="K_pneumoniae")

species_dnds_dids_intergenic_annotation_data_long_summary_wide_M_tuberculosis <- dcast(species_dnds_dids_intergenic_annotation_data_long_summary_wide, dS_bin ~ Category, value.var="M_tuberculosis")


wilcox.test(species_dnds_dids_intergenic_annotation_data_long_summary_wide_E_coli$dI.dS_promoter, species_dnds_dids_intergenic_annotation_data_long_summary_wide_E_coli$dI.dS_terminator, alternative="less")

wilcox.test(species_dnds_dids_intergenic_annotation_data_long_summary_wide_S_aureus$dI.dS_promoter, species_dnds_dids_intergenic_annotation_data_long_summary_wide_S_aureus$dI.dS_terminator, alternative="less")

wilcox.test(species_dnds_dids_intergenic_annotation_data_long_summary_wide_S_enterica$dI.dS_promoter, species_dnds_dids_intergenic_annotation_data_long_summary_wide_S_enterica$dI.dS_terminator, alternative="less")

wilcox.test(species_dnds_dids_intergenic_annotation_data_long_summary_wide_S_pneumoniae$dI.dS_promoter, species_dnds_dids_intergenic_annotation_data_long_summary_wide_S_pneumoniae$dI.dS_terminator, alternative="less")

wilcox.test(species_dnds_dids_intergenic_annotation_data_long_summary_wide_K_pneumoniae$dI.dS_promoter, species_dnds_dids_intergenic_annotation_data_long_summary_wide_K_pneumoniae$dI.dS_terminator, alternative="less")

wilcox.test(species_dnds_dids_intergenic_annotation_data_long_summary_wide_M_tuberculosis$dI.dS_promoter, species_dnds_dids_intergenic_annotation_data_long_summary_wide_M_tuberculosis$dI.dS_terminator, alternative="less")

facet_labels=c(E_coli="E. coli", S_aureus="S. aureus", S_enterica="S. enterica", S_pneumoniae="S. pneumoniae", K_pneumoniae="K. pneumoniae", M_tuberculosis="M. tuberculosis")
category_breaks=c("dN.dS", "dI.dS_rbs", "dI.dS_non_coding_RNA", "dI.dS_promoter", "dI.dS_terminator", "dI.dS_unannotated")
category_labels=c("dN/dS", "dI/dS\nRBS", "dI/dS\nNon coding RNA", "dI/dS\nPromoter", "dI/dS\nTerminator", "dI/dS\nUnannotated")

pairwise_dnds_dids_intergenic_annotation_plot <- ggplot() +
  geom_boxplot(data=species_dnds_dids_intergenic_annotation_data_long_summary, aes(x=Category, y=dX.dS), outlier.size=NA) +
  geom_point(data=species_dnds_dids_intergenic_annotation_data_long_summary, aes(x=Category, y=dX.dS, colour=Category), position=position_jitter(w=0.6)) +
  coord_cartesian(ylim=c(0, 4)) +
  facet_wrap(~Species, ncol=3, labeller=labeller(Species=facet_labels)) +
  scale_x_discrete(breaks=category_breaks, labels=category_labels) +
  labs(y="dX/dS") +
  theme(legend.position="none",
        strip.text.x=element_text(face="italic"))


#species_dnds_dids_intergenic_annotation_data_long_CC <- species_dnds_dids_intergenic_annotation_data_long_summary[ which(species_dnds_dids_intergenic_annotation_data_long_summary$dS_bin < 0.001), ]
species_dnds_dids_intergenic_annotation_data_long_CC <- species_dnds_dids_intergenic_annotation_data_long[ which(species_dnds_dids_intergenic_annotation_data_long$dS < 0.001), ]

#ggplot(species_dnds_dids_intergenic_annotation_data_long_CC, aes(x=Category, y=dX.dS, colour=Category)) +
#  geom_boxplot() +
#  scale_y_continuous(limits=c(0, 2)) +
#  facet_wrap(~Species, ncol=2)

#out_file_pdf <- paste(base_dir, "/Figures/Figure_5", ".pdf", sep="")
out_file_tif <- paste(base_dir, "/Figures/Figure_5", ".tif", sep="")

#pdf(file=out_file_pdf, height=10, width=15)
tiff(file=out_file_tif, height=10, width=15, units="in", res=100)

pairwise_dnds_dids_intergenic_annotation_plot

dev.off()

#####

species_dnds_dids_intergenic_annotation_data_long_2 <- melt(species_dnds_dids_intergenic_annotation_data, measure.vars=c("dI.dS_promoter", "dI.dS_terminator", "dI.dS_non_coding_RNA", "dI.dS_unannotated"), variable.name="Category", value.name="dX.dS")

species_dnds_dids_intergenic_annotation_data_long_2$dS_bin <- cut(species_dnds_dids_intergenic_annotation_data_long_2$dS, breaks=seq(0, 0.1, 0.0001), labels=seq(0.0001, 0.1, 0.0001))

species_dnds_dids_intergenic_annotation_data_long_2_summary <- summarySE(species_dnds_dids_intergenic_annotation_data_long_2, measurevar="dX.dS", groupvars=c("Species", "Category"))

species_dnds_dids_intergenic_annotation_data_long_2_summary$Category_2 <- NA

row_count <- nrow(species_dnds_dids_intergenic_annotation_data_long_2_summary)

for(i in 1:row_count){
  if(species_dnds_dids_intergenic_annotation_data_long_2_summary$Category[i] == "dI.dS_promoter"){
    species_dnds_dids_intergenic_annotation_data_long_2_summary$Category_2[i] <- "Promoter"
  }else if(species_dnds_dids_intergenic_annotation_data_long_2_summary$Category[i] == "dI.dS_terminator"){
    species_dnds_dids_intergenic_annotation_data_long_2_summary$Category_2[i] <- "Terminator"
  }else if(species_dnds_dids_intergenic_annotation_data_long_2_summary$Category[i] == "dI.dS_non_coding_RNA"){
    species_dnds_dids_intergenic_annotation_data_long_2_summary$Category_2[i] <- "Non_coding_RNA"
  }else if(species_dnds_dids_intergenic_annotation_data_long_2_summary$Category[i] == "dI.dS_unannotated"){
    species_dnds_dids_intergenic_annotation_data_long_2_summary$Category_2[i] <- "Unannotated"
  }
}

species_dnds_dids_intergenic_annotation_data_long_2_summary <- arrange(species_dnds_dids_intergenic_annotation_data_long_2_summary, Species, Category_2)

#####

file=paste(base_dir, file_base_1, "M_tuberculosis", file_base_2, "M_tuberculosis", file_base_3, sep="")

M_tuberculosis_dnds_dids_intergenic_annotation_data <- read.csv(file=file, header=TRUE)

M_tuberculosis_promoter_plot <- ggplot(M_tuberculosis_dnds_dids_intergenic_annotation_data, aes(x=dI.dS_promoter)) +
  geom_histogram(binwidth=0.3) +
  scale_x_continuous(breaks=seq(0, 8, 1)) +
  labs(x="dI/dS", y="Count")

#out_file_pdf <- paste(base_dir, "/Figures/Figure_S6", ".pdf", sep="")
out_file_tif <- paste(base_dir, "/Figures/Figure_S6", ".tif", sep="")

#pdf(file=out_file_pdf, height=10, width=15)
tiff(file=out_file_tif, height=10, width=15, units="in", res=100)

M_tuberculosis_promoter_plot

dev.off()
