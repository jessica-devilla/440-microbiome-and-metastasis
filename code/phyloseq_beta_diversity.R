
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("phyloseq")
#https://joey711.github.io/phyloseq/import-data.html

#BiocManager::install("mia")
#https://microbiome.github.io/mia/articles/mia.html

#BiocManager::install("microbiome")


suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(phyloseq)
  library(ggplot2)
  library(mia)
  library(scater)
  library(umap)
  library(microbiome)
  library("plyr")
})


physeq_beta_diversity <- function(physeq, dist_methods, name){
  
  physeq <- subset_taxa(physeq, Genus != "-1")
  plist <- vector("list", length(dist_methods))
  names(plist) = dist_methods
  
  distance_matrices <- list()
  
  for( i in dist_methods ){
    # Calculate distance matrix
    iDist <- phyloseq::distance(physeq, method=i)
    # Calculate ordination
    iMDS  <- ordinate(physeq, "MDS", distance=iDist)
    ## Make plot
    # Don't carry over previous plot (if error, p will be blank)
    p <- NULL
    # Create plot, store as temp variable, p
    p <- plot_ordination(physeq, iMDS, color="pathologic_stage_label", shape="data_submitting_center_label")
    # Add title to each plot
    #p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
    
    # Save the graphic to file.
    plist[[i]] = p
    
    # Store the distance matrix
    distance_matrices[[i]] <- iDist
  }
  
  
  df = ldply(plist, function(x) x$data)
  names(df)[1] <- "distance"
  p = ggplot(df, aes(Axis.1, Axis.2, color=pathologic_stage_label, shape=data_submitting_center_label))+ 
    geom_point(size=2, alpha=0.75) + theme_minimal() +
    facet_wrap(~distance, scales="free") +
    ggtitle("MDS on Bray-Curtus Distance") +
    xlab("NMDS1") +ylab("NMDS2") +
    scale_color_manual(values = c("#CED5F3", "#8070FE", "#EAB606","#FC4703"))+
    labs(color="Stage", shape ="Submitting Center")+
    theme(legend.title = element_text(size = 10),  # Set legend title size
          legend.text = element_text(size = 8))
  
  #print(p)
  filename <- paste0("figures/beta_diversity/phyloseq_beta_diversity_", name, ".png")
  ggsave(filename, plot = p,width=7, height=5)
  
  return(distance_matrices)
}


## Sample run code
#physeq <- readRDS('data/coad_raw_UNC_phyloseq.rds')
#dist_methods <- unlist(distanceMethodList)
#print(dist_methods)
#dist_methods <- c("bray")
#distance_matrixes <- physeq_beta_diversity(physeq, dist_methods,"raw_UNC_bray")

#bray <- distance_matrixes$bray


