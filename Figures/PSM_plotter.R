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

file_base_1="/Analysis/Mutation/"
file_base_2="_Mutation/"
file_base_3="_PSM.csv"

threshold_array <- c("threshold_0/", "", "threshold_99/")
threshold_breaks <- c("relaxed_core", "core", "strict_core")
threshold_labels <- c("Relaxed core", "Core", "Strict core")

threshold_count <- length(threshold_array)

species_PSM_data <- NULL

for(i in 1:species_count){
  species=species_array[i]
  
  for(j in 1:threshold_count){
    
    file=paste(base_dir, file_base_1, species, file_base_2, threshold_array[j], species, file_base_3, sep="")
    
    PSM_data <- read.csv(file=file, header=TRUE)
    
    rows <- nrow(PSM_data)
    
    Species <- rep(species, rows)
    
    PSM_data <- cbind(Species, PSM_data)
    
    PSM_data$Non_singletons <- PSM_data$Total_SNPs - PSM_data$Singletons
    
    PSM_data_wide_Singletons <- dcast(PSM_data, Species ~ Category, value.var="Singletons")
    
    PSM_data_wide_Non_singletons <- dcast(PSM_data, Species ~ Category, value.var="Non_singletons")
    
    S <- c(rep("S_INTERGENIC", PSM_data_wide_Singletons$Intergenic), rep("S_SYNONYMOUS", PSM_data_wide_Singletons$Synonymous), rep("S_NONSYNONYMOUS", PSM_data_wide_Singletons$Nonsynonymous), rep("S_NONSENSE", PSM_data_wide_Singletons$Nonsense))
    N <- c(rep("N_INTERGENIC", PSM_data_wide_Non_singletons$Intergenic), rep("N_SYNONYMOUS", PSM_data_wide_Non_singletons$Synonymous), rep("N_NONSYNONYMOUS", PSM_data_wide_Non_singletons$Nonsynonymous), rep("N_NONSENSE", PSM_data_wide_Non_singletons$Nonsense))
    
    NS <- c(N, S)
    
    Intergenic <- NULL
    Synonymous <- NULL
    Nonsynonymous <- NULL
    Nonsense <- NULL
    
    for(rep in 1:1000){
      rep_NS <- sample(NS, size=length(NS), replace=TRUE)
      
      Intergenic[rep] <- length(grep("S_INTERGENIC", rep_NS)) / (length(grep("S_INTERGENIC", rep_NS)) + length(grep("N_INTERGENIC", rep_NS)))
      Synonymous[rep] <- length(grep("S_SYNONYMOUS", rep_NS)) / (length(grep("S_SYNONYMOUS", rep_NS)) + length(grep("N_SYNONYMOUS", rep_NS)))
      Nonsynonymous[rep] <- length(grep("S_NONSYNONYMOUS", rep_NS)) / (length(grep("S_NONSYNONYMOUS", rep_NS)) + length(grep("N_NONSYNONYMOUS", rep_NS)))
      Nonsense[rep] <- length(grep("S_NONSENSE", rep_NS)) / (length(grep("S_NONSENSE", rep_NS)) + length(grep("N_NONSENSE", rep_NS)))
    }
    
    rows <- NROW(Intergenic)
    
    Species <- rep(species, rows)
    
    Threshold <- rep(threshold_breaks[j], rows)
    
    resampled_PSM_data <- data.frame(Species, Threshold, Synonymous, Intergenic, Nonsynonymous, Nonsense)
    
    species_PSM_data <- rbind(species_PSM_data, resampled_PSM_data)
  }
}

species_PSM_data <- melt(species_PSM_data, id.vars=c("Species", "Threshold"), variable.name="Category", value.name="PSM")

species_PSM_data$Order <- species_PSM_data$Total_SNPs

row_count <- nrow(species_PSM_data)

for(i in 1:row_count){
  if(species_PSM_data$Category[i] == "Synonymous"){
    species_PSM_data$Order[i] <- 1
  }else if(species_PSM_data$Category[i] == "Intergenic"){
    species_PSM_data$Order[i] <- 2
  }else if(species_PSM_data$Category[i] == "Nonsynonymous"){
    species_PSM_data$Order[i] <- 3
  }else if(species_PSM_data$Category[i] == "Nonsense"){
    species_PSM_data$Order[i] <- 4
  }
}

species_PSM_data_summary <- summarySE(species_PSM_data, groupvars=c("Species", "Threshold", "Category", "Order"), measurevar="PSM")

facet_labels=c(E_coli="E. coli", S_aureus="S. aureus", S_enterica="S. enterica", S_pneumoniae="S. pneumoniae", K_pneumoniae="K. pneumoniae", M_tuberculosis="M. tuberculosis")
category_breaks=c("Synonymous", "Intergenic", "Nonsynonymous", "Nonsense")
category_labels=c("Synonymous", "Intergenic", "Non-synonymous", "Nonsense")

PSM_plot <- ggplot(species_PSM_data_summary, aes(x=reorder(Category, Order), y=PSM, fill=Threshold)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(x=reorder(Category, Order), ymin=PSM-ci, ymax=PSM+ci, width=0.5), position=position_dodge(width=0.9)) +
  scale_y_continuous(limits=c(0, 0.8)) +
  facet_wrap(~Species, ncol=3, labeller=labeller(Species=facet_labels)) +
  scale_x_discrete(breaks=category_breaks, labels=category_labels) +
  labs(x="Category", y="PSM") +
  scale_fill_discrete(breaks=threshold_breaks, labels=threshold_labels) +
  theme(strip.text.x=element_text(face="italic"),
        axis.text.x=element_text(size=10))

#out_file_pdf <- paste(base_dir, "/Figures/Figure_2", ".pdf", sep="")
out_file_tif <- paste(base_dir, "/Figures/Figure_2", ".tif", sep="")

#pdf(file=out_file_pdf, height=10, width=15)
tiff(file=out_file_tif, height=10, width=15, units="in", res=100)

PSM_plot

dev.off()
