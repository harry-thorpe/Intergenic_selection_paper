library(cowplot)
library(reshape2)

species_array=c("E_coli", "S_aureus", "S_enterica", "S_pneumoniae", "K_pneumoniae", "M_tuberculosis")
species_count=length(species_array)

file_base_1="/media/harry/extra/Intergenic_variation_paper/Analysis/Sequence_summary/"
file_base_2="_Sequence_summary/"
file_base_3="_site_counts.csv"

species_site_data <- NULL

for(i in 1:species_count){
 # if(i == 1){
    species=species_array[i]
    
    file=paste(file_base_1, species, file_base_2, species, file_base_3, sep="")
    
    site_data <- read.csv(file=file, header=TRUE)
    
    site_data <- site_data[site_data$Category != "Nonsense", ]
    
    site_data$Total_GC_content <- ((sum(site_data$G) + sum(site_data$C)) / sum(site_data$Total))
    
    rows <- nrow(site_data)
    
    Species <- rep(species, rows)
    
    site_data <- cbind(Species, site_data)
    
    species_site_data <- rbind(species_site_data, site_data)
#  }
}

#species_site_data <- species_site_data[species_site_data$Category != "Nonsense", ]

species_site_plot <- ggplot(species_site_data, aes(x=Total_GC_content, y=GC_content, colour=Category)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  scale_x_continuous(limits=c(0.3,0.7), breaks=seq(0.3, 0.7, 0.1)) +
  scale_y_continuous(limits=c(0,0.9), breaks=seq(0, 0.9, 0.1)) +
  coord_equal(ratio=0.5) +
  labs(x="Genome GC content", y="Category GC content")

#pdf(file="/media/harry/extra/Intergenic_variation_paper/Figures/Species_site_GC.pdf", height=10, width=15)
tiff(file="/media/harry/extra/Intergenic_variation_paper/Figures/Species_site_GC.tif", height=10, width=15, units="in", res=100)

species_site_plot

dev.off()
