
rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(phyloseq); packageVersion("phyloseq")
  library(ggplot2)
  library(mia)
  library(scater)
  library(umap)
  library(microbiome)
  library(plyr)
})

source("code/clean_kraken_data.R")
source("code/run_norm_func.R")
source("code/make_phyloseq_obj.R")
source("code/phyloseq_beta_diversity.R")


# import and process datasets of interest

#RAW DATA - NORMS NOT WORKING YET
kraken_meta <- readRDS("data/kraken_metaCOAD.RDS")
kraken_data <- readRDS("data/kraken_COAD_raw.RDS")

# UNC DATA

#kraken_meta <- readRDS("data/kraken_meta_norm_filtered.RDS")
#kraken_data<- readRDS("data/kraken_norm_filtered.RDS") # norm'd
#kraken_data <-readRDS("data/kraken_raw_filtered.RDS") # raw

result <- clean_kraken_data(kraken_data, kraken_meta)
kraken_data <- result$kraken_data
kraken_meta <- result$kraken_meta

kraken_data_t <- t(kraken_data)

## first run the normalization methods on each dataset
# List of normalization methods
norm_methods <- c("DeSEQ", "RLE+", "RLE_poscounts", "TSS", "UQ", "CSS",'CLR_poscounts', "logcpm", "rarefy", "CLR+", "MED", "GMPR", "CLR_poscounts")


# Initialize a list to store normalized data frames
normalized_dataframes <- list()

# Loop through each normalization method
for (method in norm_methods) {
  # Try applying normalization function to kraken_data
  tryCatch({
    # Print error message if normalization fails for a method
    cat(paste("Running method: ", method, "\n"))
    # Apply normalization function
    normalized_data <- norm.func(kraken_data_t, method)
    normalized_data <- as.data.frame(normalized_data)
    
    
    # Store the transposed normalized dataframe in the list with method as the key
    normalized_dataframes[[method]] <- normalized_data
  }, error = function(e) {
    # Print error message if normalization fails for a method
    cat(paste("Error occurred for method", method, ":", conditionMessage(e), "\n"))
  })
}

## then iterate through each normalized dataset
# create phyloseq object using data and metadata

# Initialize a list to store phyloseq objects
phyloseq_objects <- list()

# Loop through each normalized data frame
for (method in names(normalized_dataframes)) {
  cat(paste("Running method: ", method, "\n"))
  # Create phyloseq object using create_phyloseq function
  physeq <- make_phyloseq_object(normalized_dataframes[[method]], kraken_meta)
  # Store the phyloseq object in the list with method as the key
  phyloseq_objects[[method]] <- physeq
  
  filename <- paste0("allsamples_", method)
  # run beta diversity function and make plot
  distance_matrixes <- physeq_beta_diversity(physeq, dist_methods = c("bray"), name = filename)
  print(mean(bray <- distance_matrixes$bray))
  
}

# plot zico 1-4 for all normalized
# calculate beta diversity for all normalized and plot
# calculate alpha diversity for all normalized and plot

# calculate some metric of similarity between the normalized and the raw matrices and show as heat map

