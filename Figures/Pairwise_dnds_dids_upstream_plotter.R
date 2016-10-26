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

file_base_1="/Analysis/Pairwise_dnds_dids_upstream/"
file_base_2="_Pairwise_dnds_dids_upstream/"
file_base_3="_dnds_dids_upstream.csv"

file_base_1_1="/Analysis/Pairwise_dnds_dids_upstream_mutation_bias_correction/"
file_base_2_2="_Pairwise_dnds_dids_upstream_mutation_bias_correction/"
file_base_3_3="_dnds_dids_intergenic_upstream_simulated.csv"

species_dnds_dids_upstream_data <- NULL

for(i in 1:species_count){
  species=species_array[i]
  
  file=paste(base_dir, file_base_1_1, species, file_base_2_2, species, file_base_3_3, sep="")
  
  dnds_dids_upstream_data <- read.csv(file=file, header=TRUE)
  
  mean_dN.dS <- mean(dnds_dids_upstream_data$dN.dS)
  mean_dI.dS <- mean(dnds_dids_upstream_data$dI.dS)
  mean_dI.dS_upstream_30 <- mean(dnds_dids_upstream_data$dI.dS_upstream_30)
  
  file=paste(base_dir, file_base_1, species, file_base_2, species, file_base_3, sep="")
  
  dnds_dids_upstream_data <- read.csv(file=file, header=TRUE)
  
  dnds_dids_upstream_data$dN.dS <- dnds_dids_upstream_data$dN.dS/mean_dN.dS
  dnds_dids_upstream_data$dI.dS <- dnds_dids_upstream_data$dI.dS/mean_dI.dS
  dnds_dids_upstream_data$dI.dS_upstream_30 <- dnds_dids_upstream_data$dI.dS_upstream_30/mean_dI.dS_upstream_30
  
  rows <- nrow(dnds_dids_upstream_data)
  
  Species <- rep(species, rows)
  
  dnds_dids_upstream_data <- cbind(Species, dnds_dids_upstream_data)
  
  print(wilcox.test(dnds_dids_upstream_data$dI.dS_upstream_30, dnds_dids_upstream_data$dI.dS, alternative="less"))
  
  species_dnds_dids_upstream_data <- rbind(species_dnds_dids_upstream_data, dnds_dids_upstream_data)
}

species_dnds_dids_upstream_data_long <- melt(species_dnds_dids_upstream_data, measure.vars=c("dN.dS", "dI.dS_upstream_30", "dI.dS"), variable.name="Category", value.name="dX.dS")

species_dnds_dids_upstream_data_long$dS_bin <- cut(species_dnds_dids_upstream_data_long$dS, breaks=seq(0, 0.1, 0.0001), labels=seq(0.0001, 0.1, 0.0001))

species_dnds_dids_upstream_data_long_summary <- summarySE(species_dnds_dids_upstream_data_long, measurevar="dX.dS", groupvars=c("Species", "Category", "dS_bin"))

facet_labels=c(E_coli="E. coli", S_aureus="S. aureus", S_enterica="S. enterica", S_pneumoniae="S. pneumoniae", K_pneumoniae="K. pneumoniae", M_tuberculosis="M. tuberculosis")
category_breaks=c("dN.dS", "dI.dS_upstream_30", "dI.dS")
category_labels=c("dN/dS", "dI/dS\nUpstream", "dI/dS\nAll")

pairwise_dnds_dids_upstream_plot <- ggplot() +
  geom_boxplot(data=species_dnds_dids_upstream_data_long_summary, aes(x=Category, y=dX.dS), outlier.size=NA) +
  geom_point(data=species_dnds_dids_upstream_data_long_summary, aes(x=Category, y=dX.dS, colour=Category), position=position_jitter(w=0.6)) +
  coord_cartesian(ylim=c(0,2)) +
  facet_wrap(~Species, ncol=3, labeller=labeller(Species=facet_labels)) +
  scale_x_discrete(breaks=category_breaks, labels=category_labels) +
  labs(y="dX/dS") +
  theme(legend.position="none",
        strip.text.x=element_text(face="italic"))

#out_file_pdf <- paste(base_dir, "/Figures/Figure_S3", ".pdf", sep="")
out_file_tif <- paste(base_dir, "/Figures/Figure_S3", ".tif", sep="")

#pdf(file=out_file_pdf, height=10, width=15)
tiff(file=out_file_tif, height=10, width=15, units="in", res=100)

pairwise_dnds_dids_upstream_plot

dev.off()
