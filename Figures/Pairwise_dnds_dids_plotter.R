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

file_base_1="/Analysis/Pairwise_dnds_dids/"
file_base_2="_Pairwise_dnds_dids/"
file_base_3="_dnds_dids.csv"

file_base_1_1="/Analysis/Pairwise_dnds_dids_mutation_bias_correction/"
file_base_2_2="_Pairwise_dnds_dids_mutation_bias_correction/"
file_base_3_3="_dnds_dids_simulated.csv"

species_dnds_dids_data <- NULL

for(i in 1:species_count){
  species=species_array[i]
  
  file=paste(base_dir, file_base_1_1, species, file_base_2_2, species, file_base_3_3, sep="")
  
  dnds_dids_data <- read.csv(file=file, header=TRUE)
  
  mean_dN.dS <- mean(dnds_dids_data$dN.dS)
  mean_dI.dS <- mean(dnds_dids_data$dI.dS)
  
  file=paste(base_dir, file_base_1, species, file_base_2, species, file_base_3, sep="")
  
  dnds_dids_data <- read.csv(file=file, header=TRUE)
  
  dnds_dids_data$dN.dS <- dnds_dids_data$dN.dS/mean_dN.dS
  dnds_dids_data$dI.dS <- dnds_dids_data$dI.dS/mean_dI.dS
  
  rows <- nrow(dnds_dids_data)
  
  print(cor.test(dnds_dids_data$dN.dS, dnds_dids_data$dS, method="spearman"))
  print(cor.test(dnds_dids_data$dI.dS, dnds_dids_data$dS, method="spearman"))
  
  Species <- rep(species, rows)
  
  dnds_dids_data <- cbind(Species, dnds_dids_data)
  
  species_dnds_dids_data <- rbind(species_dnds_dids_data, dnds_dids_data)
}

species_dnds_dids_data_long <- melt(species_dnds_dids_data, measure.vars=c("dN.dS", "dI.dS"), variable.name="Category", value.name="dX.dS")

species_dnds_dids_data_long$dS_bin <- cut(species_dnds_dids_data_long$dS, breaks=seq(0, 0.1, 0.0001), labels=seq(0.0001, 0.1, 0.0001))

species_dnds_dids_data_long$dS_bin <- as.numeric(as.character(species_dnds_dids_data_long$dS_bin))

species_dnds_dids_data_long_summary <- summarySE(species_dnds_dids_data_long, measurevar="dX.dS", groupvars=c("Species", "Category", "dS_bin"))

pairwise_dnds_dids_summary_plot <- ggplot(species_dnds_dids_data_long_summary, aes(x=dS_bin, y=dX.dS, colour=Category)) +
  geom_point() +
  geom_errorbar(aes(ymin=dX.dS-se, ymax=dX.dS+se)) +
  coord_cartesian(ylim=c(0,2)) +
  facet_wrap(~Species, ncol=2, scales="free")


species_dnds_dids_data_long_summary_wide <- dcast(species_dnds_dids_data_long_summary, Species + dS_bin ~ Category, value.var="dX.dS")

species_dnds_dids_data_long_summary_wide_dN.dS <- dcast(species_dnds_dids_data_long_summary_wide, dS_bin ~ Species, value.var="dN.dS")

species_dnds_dids_data_long_summary_wide_dI.dS <- dcast(species_dnds_dids_data_long_summary_wide, dS_bin ~ Species, value.var="dI.dS")

cor.test(species_dnds_dids_data_long_summary_wide_dN.dS$E_coli, species_dnds_dids_data_long_summary_wide_dN.dS$dS_bin, method="spearman")
cor.test(species_dnds_dids_data_long_summary_wide_dI.dS$E_coli, species_dnds_dids_data_long_summary_wide_dI.dS$dS_bin, method="spearman")

cor.test(species_dnds_dids_data_long_summary_wide_dN.dS$S_aureus, species_dnds_dids_data_long_summary_wide_dN.dS$dS_bin, method="spearman")
cor.test(species_dnds_dids_data_long_summary_wide_dI.dS$S_aureus, species_dnds_dids_data_long_summary_wide_dI.dS$dS_bin, method="spearman")

cor.test(species_dnds_dids_data_long_summary_wide_dN.dS$S_enterica, species_dnds_dids_data_long_summary_wide_dN.dS$dS_bin, method="spearman")
cor.test(species_dnds_dids_data_long_summary_wide_dI.dS$S_enterica, species_dnds_dids_data_long_summary_wide_dI.dS$dS_bin, method="spearman")

cor.test(species_dnds_dids_data_long_summary_wide_dN.dS$S_pneumoniae, species_dnds_dids_data_long_summary_wide_dN.dS$dS_bin, method="spearman")
cor.test(species_dnds_dids_data_long_summary_wide_dI.dS$S_pneumoniae, species_dnds_dids_data_long_summary_wide_dI.dS$dS_bin, method="spearman")

cor.test(species_dnds_dids_data_long_summary_wide_dN.dS$K_pneumoniae, species_dnds_dids_data_long_summary_wide_dN.dS$dS_bin, method="spearman")
cor.test(species_dnds_dids_data_long_summary_wide_dI.dS$K_pneumoniae, species_dnds_dids_data_long_summary_wide_dI.dS$dS_bin, method="spearman")

cor.test(species_dnds_dids_data_long_summary_wide_dN.dS$M_tuberculosis, species_dnds_dids_data_long_summary_wide_dN.dS$dS_bin, method="spearman")
cor.test(species_dnds_dids_data_long_summary_wide_dI.dS$M_tuberculosis, species_dnds_dids_data_long_summary_wide_dI.dS$dS_bin, method="spearman")

facet_labels=c(E_coli="E. coli", S_aureus="S. aureus", S_enterica="S. enterica", S_pneumoniae="S. pneumoniae", K_pneumoniae="K. pneumoniae", M_tuberculosis="M. tuberculosis")
category_breaks=c("dN.dS", "dI.dS")
category_labels=c("dN/dS", "dI/dS")

pairwise_dnds_dids_plot <- ggplot(species_dnds_dids_data_long, aes(x=dS, y=dX.dS, colour=Category)) +
  geom_point() +
  coord_cartesian(ylim=c(0, 2)) +
  facet_wrap(~Species, ncol=3, scales="free", labeller=labeller(Species=facet_labels)) +
  scale_colour_discrete(breaks=category_breaks, labels=category_labels) +
  #guides(colour=guide_legend(override.aes=list(alpha=1))) +
  labs(y="dX/dS") +
  theme(strip.text.x=element_text(face="italic"))

#out_file_pdf <- paste(base_dir, "/Figures/Figure_S2", ".pdf", sep="")
out_file_tif <- paste(base_dir, "/Figures/Figure_S2", ".tif", sep="")

#pdf(file=out_file_pdf, height=10, width=15)
tiff(file=out_file_tif, height=10, width=15, units="in", res=100)

pairwise_dnds_dids_plot

dev.off()

#####

species_dnds_dids_data_long$Comparison <- cut(species_dnds_dids_data_long$dS, breaks=c(0, 0.001, 1), labels=c("Within_clonal_complex", "Between_clonal_complex"))

facet_labels=c(E_coli="E. coli", S_aureus="S. aureus", S_enterica="S. enterica", S_pneumoniae="S. pneumoniae", K_pneumoniae="K. pneumoniae", M_tuberculosis="M. tuberculosis")
category_breaks=c("dN.dS", "dI.dS")
category_labels=c("dN/dS", "dI/dS")
comparison_breaks=c("Within_clonal_complex", "Between_clonal_complex")
comparison_labels=c("Within CC", "Between CC")

pairwise_dnds_dids_clonal_plot <- ggplot() +
  geom_boxplot(data=species_dnds_dids_data_long, aes(x=Category, y=dX.dS, colour=Comparison), notch=TRUE) +
  coord_cartesian(ylim=c(0,2)) +
  facet_grid(.~Species, labeller=labeller(Species=facet_labels)) +
  scale_colour_discrete(breaks=comparison_breaks, labels=comparison_labels) +
  #scale_x_discrete(breaks=category_breaks, labels=category_labels) +
  labs(y="dX/dS") +
  theme(strip.text.x=element_text(face="italic"))

#out_file_pdf <- paste(base_dir, "/Figures/Figure_3", ".pdf", sep="")
out_file_tif <- paste(base_dir, "/Figures/Figure_3", ".tif", sep="")

#pdf(file=out_file_pdf, height=10, width=15)
tiff(file=out_file_tif, height=10, width=15, units="in", res=100)

pairwise_dnds_dids_clonal_plot

dev.off()

species_dnds_dids_data_wide <- dcast(species_dnds_dids_data_long, Species + Isolate_1 + Isolate_2 + Comparison + dS + dS_bin ~ Category, value.var="dX.dS")

species_dnds_dids_data_wide_dN <- dcast(species_dnds_dids_data_wide, Species + Isolate_1 + Isolate_2 + dS + dS_bin ~ Comparison, value.var="dN.dS")

species_dnds_dids_data_wide_dI <- dcast(species_dnds_dids_data_wide, Species + Isolate_1 + Isolate_2 + dS + dS_bin ~ Comparison, value.var="dI.dS")

species_dnds_dids_data_wide_dN_within <- dcast(species_dnds_dids_data_wide_dN, Isolate_1 + Isolate_2 + dS + dS_bin ~ Species, value.var="Within_clonal_complex")

species_dnds_dids_data_wide_dN_between <- dcast(species_dnds_dids_data_wide_dN, Isolate_1 + Isolate_2 + dS + dS_bin ~ Species, value.var="Between_clonal_complex")

species_dnds_dids_data_wide_dI_within <- dcast(species_dnds_dids_data_wide_dI, Isolate_1 + Isolate_2 + dS + dS_bin ~ Species, value.var="Within_clonal_complex")

species_dnds_dids_data_wide_dI_between <- dcast(species_dnds_dids_data_wide_dI, Isolate_1 + Isolate_2 + dS + dS_bin ~ Species, value.var="Between_clonal_complex")

wilcox.test(species_dnds_dids_data_wide_dN_within$E_coli, species_dnds_dids_data_wide_dN_between$E_coli, alternative="greater")
wilcox.test(species_dnds_dids_data_wide_dI_within$E_coli, species_dnds_dids_data_wide_dI_between$E_coli, alternative="greater")

wilcox.test(species_dnds_dids_data_wide_dI_within$E_coli, mu=1, alternative="less")

wilcox.test(species_dnds_dids_data_wide_dN_within$S_aureus, species_dnds_dids_data_wide_dN_between$S_aureus, alternative="greater")
wilcox.test(species_dnds_dids_data_wide_dI_within$S_aureus, species_dnds_dids_data_wide_dI_between$S_aureus, alternative="greater")

wilcox.test(species_dnds_dids_data_wide_dI_within$S_aureus, mu=1, alternative="less")

wilcox.test(species_dnds_dids_data_wide_dN_within$S_enterica, species_dnds_dids_data_wide_dN_between$S_enterica, alternative="greater")
wilcox.test(species_dnds_dids_data_wide_dI_within$S_enterica, species_dnds_dids_data_wide_dI_between$S_enterica, alternative="greater")

wilcox.test(species_dnds_dids_data_wide_dI_within$S_enterica, mu=1, alternative="less")

wilcox.test(species_dnds_dids_data_wide_dN_within$S_pneumoniae, species_dnds_dids_data_wide_dN_between$S_pneumoniae, alternative="greater")
wilcox.test(species_dnds_dids_data_wide_dI_within$S_pneumoniae, species_dnds_dids_data_wide_dI_between$S_pneumoniae, alternative="greater")

wilcox.test(species_dnds_dids_data_wide_dI_within$S_pneumoniae, mu=1, alternative="less")

wilcox.test(species_dnds_dids_data_wide_dN_within$K_pneumoniae, species_dnds_dids_data_wide_dN_between$K_pneumoniae, alternative="greater")
wilcox.test(species_dnds_dids_data_wide_dI_within$K_pneumoniae, species_dnds_dids_data_wide_dI_between$K_pneumoniae, alternative="greater")

wilcox.test(species_dnds_dids_data_wide_dI_within$K_pneumoniae, mu=1, alternative="less")
