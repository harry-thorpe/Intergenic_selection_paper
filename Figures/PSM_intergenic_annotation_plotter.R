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

file_base_1="/Analysis/Mutation_intergenic_annotation/"
file_base_2="_Mutation_intergenic_annotation/"
file_base_3="_PSM_intergenic.csv"

species_PSM_intergenic_annotation_data <- NULL

for(i in 1:species_count){
  species=species_array[i]
  
  file=paste(base_dir, file_base_1, species, file_base_2, species, file_base_3, sep="")
  
  PSM_intergenic_annotation_data <- read.csv(file=file, header=TRUE)
  
  rows <- nrow(PSM_intergenic_annotation_data)
  
  Species <- rep(species, rows)
  
  PSM_intergenic_annotation_data <- cbind(Species, PSM_intergenic_annotation_data)
  
  PSM_intergenic_annotation_data$Non_singletons <- PSM_intergenic_annotation_data$Total_SNPs - PSM_intergenic_annotation_data$Singletons
  
  PSM_intergenic_annotation_data_wide_Singletons <- dcast(PSM_intergenic_annotation_data, Species ~ Category, value.var="Singletons")
  
  PSM_intergenic_annotation_data_wide_Non_singletons <- dcast(PSM_intergenic_annotation_data, Species ~ Category, value.var="Non_singletons")
  
  S <- c(rep("S_rbs", PSM_intergenic_annotation_data_wide_Singletons$rbs), rep("S_Unannotated", PSM_intergenic_annotation_data_wide_Singletons$Unannotated), rep("S_Promoter", PSM_intergenic_annotation_data_wide_Singletons$Promoter), rep("S_Terminator", PSM_intergenic_annotation_data_wide_Singletons$Terminator), rep("S_Non_coding_RNA", PSM_intergenic_annotation_data_wide_Singletons$Non_coding_RNA))
  N <- c(rep("N_rbs", PSM_intergenic_annotation_data_wide_Non_singletons$rbs), rep("N_Unannotated", PSM_intergenic_annotation_data_wide_Non_singletons$Unannotated), rep("N_Promoter", PSM_intergenic_annotation_data_wide_Non_singletons$Promoter), rep("N_Terminator", PSM_intergenic_annotation_data_wide_Non_singletons$Terminator), rep("N_Non_coding_RNA", PSM_intergenic_annotation_data_wide_Non_singletons$Non_coding_RNA))
  
  NS <- c(N, S)
  
  rbs <- NULL
  Unannotated <- NULL
  Promoter <- NULL
  Terminator <- NULL
  Non_coding_RNA <- NULL
  
  for(rep in 1:10){
    rep_NS <- sample(NS, size=length(NS), replace=TRUE)
    
    rbs[rep] <- length(grep("S_rbs", rep_NS)) / (length(grep("S_rbs", rep_NS)) + length(grep("N_rbs", rep_NS)))
    Unannotated[rep] <- length(grep("S_Unannotated", rep_NS)) / (length(grep("S_Unannotated", rep_NS)) + length(grep("N_Unannotated", rep_NS)))
    Promoter[rep] <- length(grep("S_Promoter", rep_NS)) / (length(grep("S_Promoter", rep_NS)) + length(grep("N_Promoter", rep_NS)))
    Terminator[rep] <- length(grep("S_Terminator", rep_NS)) / (length(grep("S_Terminator", rep_NS)) + length(grep("N_Terminator", rep_NS)))
    Non_coding_RNA[rep] <- length(grep("S_Non_coding_RNA", rep_NS)) / (length(grep("S_Non_coding_RNA", rep_NS)) + length(grep("N_Non_coding_RNA", rep_NS)))
  }
  
  rows <- NROW(Unannotated)
  
  Species <- rep(species, rows)
  
  resampled_PSM_intergenic_annotation_data <- data.frame(Species, Unannotated, Promoter, Terminator, Non_coding_RNA, rbs)
  
  species_PSM_intergenic_annotation_data <- rbind(species_PSM_intergenic_annotation_data, resampled_PSM_intergenic_annotation_data)
}

species_PSM_intergenic_annotation_data <- melt(species_PSM_intergenic_annotation_data, id.vars=c("Species"), variable.name="Category", value.name="PSM")

species_PSM_intergenic_annotation_data$Order <- species_PSM_intergenic_annotation_data$PSM

row_count <- nrow(species_PSM_intergenic_annotation_data)

for(i in 1:row_count){
  if(species_PSM_intergenic_annotation_data$Category[i] == "Unannotated"){
    species_PSM_intergenic_annotation_data$Order[i] <- 1
  }else if(species_PSM_intergenic_annotation_data$Category[i] == "Terminator"){
    species_PSM_intergenic_annotation_data$Order[i] <- 2
  }else if(species_PSM_intergenic_annotation_data$Category[i] == "Promoter"){
    species_PSM_intergenic_annotation_data$Order[i] <- 3
  }else if(species_PSM_intergenic_annotation_data$Category[i] == "Non_coding_RNA"){
    species_PSM_intergenic_annotation_data$Order[i] <- 4
  }else if(species_PSM_intergenic_annotation_data$Category[i] == "rbs"){
    species_PSM_intergenic_annotation_data$Order[i] <- 5
  }
}

species_PSM_intergenic_annotation_data_summary <- summarySE(species_PSM_intergenic_annotation_data, groupvars=c("Species", "Category", "Order"), measurevar="PSM")

facet_labels=c(E_coli="E. coli", S_aureus="S. aureus", S_enterica="S. enterica", S_pneumoniae="S. pneumoniae", K_pneumoniae="K. pneumoniae", M_tuberculosis="M. tuberculosis")
category_breaks=c("rbs", "Non_coding_RNA", "Promoter", "Terminator", "Unannotated")
category_labels=c("RBS", "Non coding\nRNA", "Promoter", "Terminator", "Unannotated")

PSM_intergenic_annotation_plot <- ggplot(species_PSM_intergenic_annotation_data_summary, aes(x=reorder(Category, Order), y=PSM, fill=Category)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(x=reorder(Category, Order), ymin=PSM-ci, ymax=PSM+ci), width=0.5) +
  scale_y_continuous(limits=c(0, 0.8)) +
  facet_wrap(~Species, ncol=3, labeller=labeller(Species=facet_labels)) +
  scale_x_discrete(breaks=category_breaks, labels=category_labels) +
  labs(x="Category", y="PSM") +
  theme(legend.position="none",
        strip.text.x=element_text(face="italic"),
        axis.text.x=element_text(size=10))

#out_file_pdf <- paste(base_dir, "/Figures/Figure_S6", ".pdf", sep="")
out_file_tif <- paste(base_dir, "/Figures/Figure_S6", ".tif", sep="")

#pdf(file=out_file_pdf, height=10, width=15)
tiff(file=out_file_tif, height=10, width=15, units="in", res=100)

PSM_intergenic_annotation_plot

dev.off()
