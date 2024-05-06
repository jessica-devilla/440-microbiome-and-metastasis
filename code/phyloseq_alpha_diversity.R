
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
#https://joey711.github.io/phyloseq/import-data.html

BiocManager::install("mia")
#https://microbiome.github.io/mia/articles/mia.html

BiocManager::install("microbiome")


suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(phyloseq); packageVersion("phyloseq")
  library(ggplot2)
  library(mia)
  library(scater)
  library(umap)
  library(microbiome)
})

rm(list = ls(all.names = TRUE))

source("code/make_phyloseq_obj.R")
source("code/clean_kraken_data.R")


physeq_alpha_diversity <- function(physeq){
  # NOTE - CANNOT CALCULATE ALPHA DIVERISTY WITHOUT RAW COUNTS
  alpha_div <- estimate_richness(physeq, split = TRUE, measures = c("Shannon")) #split=true calculates per sample
  p <- plot_richness(physeq, x="pathologic_stage_label", color="pathologic_stage_label",measures=c( "Shannon"))+
    labs(x = "Pathologic Stage", y = "Shannon Diversity") +
    theme(legend.position = "none")
  print(p)
  ggsave("figures/phyloseq_shannon_raw_filt.png", plot = p)
  
  # Calculate means and standard deviations of Shannon diversity for each pathologic stage level
  mean_sd_alpha_div <- aggregate(alpha_div, by = list(pathologic_stage_label = sample_data(physeq)$pathologic_stage_label), 
                                 FUN = function(x) c(mean = mean(x), sd = sd(x)))
  
  
  # Calculate upper and lower limits for error bars
  mean_sd_alpha_div$upper <- mean_sd_alpha_div$Shannon[,"mean"] + mean_sd_alpha_div$Shannon[,"sd"]
  mean_sd_alpha_div$lower <- mean_sd_alpha_div$Shannon[,"mean"] - mean_sd_alpha_div$Shannon[,"sd"]
  
  # Plot mean and standard deviation
  p2 <- ggplot(mean_sd_alpha_div, aes(x = pathologic_stage_label, y = Shannon[,"mean"]), color=patholoic_stage_label) +
    geom_bar(aes(x = pathologic_stage_label, y = Shannon[,"mean"], fill=pathologic_stage_label), stat = "identity") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.9)) +
    labs(x = "Pathologic Stage", y = "Shannon Diversity") +
    theme(legend.position = "none")
  print(p2)
  ggsave("figures/phyloseq_shannon_raw_filt_meansd.png", plot = p2)
  return(mean_sd_alpha_div)
}

mia_alpha_diversity <- function(physeq){
  # convert to tree summarized experiment using mia library and try again
  cat("making tse object \n")
  tse <- makeTreeSEFromPhyloseq(physeq)
  tse
  
  outmat <- as(otu_table(physeq), "matrix")
  taxmat <- as(tax_table(physeq), "matrix")
  
  # get most abundant feature IDs
  top_taxa <- getTopFeatures(tse, method = "mean", top = 5, assay.type = "counts")
  # Find the indices of the rows with the target OTUs in their row names
  target_rows <- which(rownames(taxmat) %in% top_taxa)
  # Print the rows
  print(taxmat[target_rows, ])
  
  # generate tidy data
  molten_data <- meltAssay(tse,
                           assay.type = "counts",
                           add_row_data = TRUE,
                           add_col_data = TRUE
  )
  
  # calculate alpha diversity metrics
  
  ## TROUBLESHOOTING WITH LOGCPM
  #index <- c("shannon","gini_simpson","inverse_simpson", "coverage", "fisher", 
  #"faith",  "log_modulo_skewness")
  # cannot calculate chao1 diveristy with non-integer values
  #cannot calculate shannon diveristy with negative values
  # fisher only accepts counts
  # faith requires row tree in input argument x
  # these 3 will run but throw 3 error about negative entries: gini-simpson, inverse-simpson, coverage
  # log-modulo-skewness is only one that works without error
  
  #Index methods
  idxs <- c("log_modulo_skewness")
  # Corresponding polished names
  names <- c(  "LogModSkewness")
  
  
  # Calculate diversities
  tse <- estimateDiversity(tse, index = idxs, name=names)
  
  # The colData contains the indices with their code names by default
  #divmat <- colData(tse)[, "LogModSkewness"]
  
  diversity_plt <- plotColData(tse, "LogModSkewness", "pathologic_stage_label",color_by = "pathologic_stage_label")+
    scale_color_manual(values = c("blue", "#8070FE", "#EAB606","#FC4703"))+
    theme(axis.text = element_text(size = 12), 
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))
  
  
  print(diversity_plt)
  
  
  
  ggsave("figures/mia_allsamples_logmodskewness_bystage.png", plot = diversity_plt, width=5.5, height=5)
  return(tse)
}

#kraken_meta <- readRDS("data/kraken_meta_norm_filtered.RDS")
#kraken_data<- readRDS("data/kraken_norm_filtered.RDS")
#kraken_data <-readRDS("data/kraken_raw_filtered.RDS")


kraken_voom_snm <- readRDS("data/kraken_COAD.RDS")
kraken_meta_voom <- readRDS("data/kraken_metaCOAD.RDS")
result_voom <- clean_kraken_data(kraken_voom_snm,kraken_meta_voom)
kraken_data <- result_voom$kraken_data
kraken_meta <- result_voom$kraken_meta


physeq <- make_phyloseq_object(kraken_data, kraken_meta)
#saveRDS(physeq, 'data/coad_raw_UNC_phyloseq.rds')

# test plot
plot_bar(physeq, x="pathologic_stage_label", fill="Phylum")


mean_sd_alpha_div <- physeq_alpha_diversity(physeq) #only works for raw counts data

tse <- mia_alpha_diversity(physeq)



###NOT USED

# Merge samples using phyloseq package
physeq_st = merge_samples(physeq, "pathologic_stage_label")
# repair variables that were damaged during merge (coerced to numeric)
sample_data(physeq_st)$pathologic_stage_label <- factor(sample_names(physeq_st))
p = plot_richness(physeq_st,x="pathologic_stage_label" ,color="pathologic_stage_label", measures=c("Shannon"))
p + geom_point(size=5, alpha=0.7)
ggsave("figures/phyloseq_shannon_raw_filt_merged.png", plot = p)



