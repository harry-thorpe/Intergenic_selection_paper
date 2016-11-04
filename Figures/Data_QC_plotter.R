library(cowplot)
library(reshape2)

# species_array=c("E_coli", "S_aureus", "S_enterica", "S_pneumoniae", "K_pneumoniae", "M_tuberculosis")
# base_dir <- "/media/harry/extra/Intergenic_selection_paper_output"

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

file_base_1="/Analysis/Data_QC/"
file_base_2="_Data_QC/"
file_base_3="_variants_subsampled.tab"

species_data_qc <- NULL

for(i in 1:species_count){
  species=species_array[i]
  
  file=paste(base_dir, file_base_1, species, file_base_2, species, file_base_3, sep="")
  
  data_qc <- read.delim(file, sep="\t", header=TRUE)
  
  rows <- nrow(data_qc)
  
  singleton_support <- nrow(data_qc[(data_qc$Support >= 0.9 & data_qc$Count == 1), ])
  non_singleton_support <- nrow(data_qc[(data_qc$Support >= 0.9 & data_qc$Count != 1), ])
  
  singleton_phred <- nrow(data_qc[(data_qc$Score >= 200 & data_qc$Count == 1), ])
  non_singleton_phred <- nrow(data_qc[(data_qc$Score >= 200 & data_qc$Count != 1), ])
  
  singleton <- nrow(data_qc[data_qc$Count == 1, ])
  non_singleton <- nrow(data_qc[data_qc$Count != 1, ])
  
  prop_singleton_support <- singleton_support/singleton
  prop_non_singleton_support <- non_singleton_support/non_singleton
  
  prop_singleton_phred <- singleton_phred/singleton
  prop_non_singleton_phred <- non_singleton_phred/non_singleton
  
  print(sprintf("%s  %f  %f  %f  %f", species, prop_singleton_support, prop_non_singleton_support, prop_singleton_phred, prop_non_singleton_phred))
  
  Species <- rep(species, rows)
  
  data_qc <- cbind(Species, data_qc)
  
  species_data_qc <- rbind(species_data_qc, data_qc)
}

# species_data_qc$Category <- "Non_singleton"
# species_data_qc$Category[species_data_qc$Count == 1] <- "Singleton"
# 
# species_breaks=c("E_coli", "S_aureus", "S_enterica", "S_pneumoniae", "K_pneumoniae", "M_tuberculosis")
# species_labels=c("E. coli", "S. aureus", "S. enterica", "S. pneumoniae", "K. pneumoniae", "M. tuberculosis")
# 
# facet_labels=c(E_coli="E. coli", S_aureus="S. aureus", S_enterica="S. enterica", S_pneumoniae="S. pneumoniae", K_pneumoniae="K. pneumoniae", M_tuberculosis="M. tuberculosis")
# 
# 
# depth_plot <- ggplot(species_data_qc, aes(x=Depth, colour=Category)) +
#   geom_density() +
#   scale_x_continuous(limits=c(0,200)) +
#   facet_wrap(~Species, ncol=6, scales="free", labeller=labeller(Species=facet_labels)) +
#   labs(x="Coverage depth", y="Density") +
#   theme(strip.text.x=element_text(face="italic"))
# 
# support_plot <- ggplot(species_data_qc, aes(x=Support, colour=Category)) +
#   geom_density() +
#   scale_x_continuous(limits=c(0.9,1.01)) +
#   facet_wrap(~Species, ncol=6, scales="free", labeller=labeller(Species=facet_labels)) +
#   labs(x="Proportion of reads supporting variant", y="Density") +
#   theme(strip.text.x=element_text(face="italic"))
# 
# phred_plot <- ggplot(species_data_qc, aes(x=Score, colour=Category)) +
#   geom_density() +
#   scale_x_continuous(limits=c(200,225)) +
#   facet_wrap(~Species, ncol=6, scales="free", labeller=labeller(Species=facet_labels)) +
#   labs(x="Phred score", y="Density") +
#   theme(strip.text.x=element_text(face="italic"))
# 
# 
# data_qc_plot <- ggdraw() +
#   draw_plot(phred_plot, 0, 0, 1, 0.33) +
#   draw_plot(support_plot, 0, 0.33, 1, 0.33) +
#   draw_plot(depth_plot, 0, 0.66, 1, 0.33) +
#   draw_plot_label(c("a", "b", "c"), c(0, 0, 0), c(1, 0.66, 0.33))
# 
# #out_file_pdf <- paste(base_dir, "/Figures/Figure_S8", ".pdf", sep="")
# out_file_tif <- paste(base_dir, "/Figures/Figure_S8", ".tif", sep="")
# 
# #pdf(file=out_file_pdf, height=7, width=15)
# tiff(file=out_file_tif, height=7, width=15, units="in", res=100)
# 
# data_qc_plot
# 
# dev.off()
