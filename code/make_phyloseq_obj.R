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


make_phyloseq_object <- function(kraken_data, kraken_meta){
  kraken_data <- as.matrix(kraken_data)
  
  
  # Extract taxonomy levels
  col_names <- colnames(kraken_data)
  num_samples <- length(col_names)
  
  split_values = strsplit(col_names, split = ".",fixed=TRUE)
  
  # Get unique prefixes
  prefixes <- unique(unlist(lapply(split_values, function(x) gsub("__.*", "", x))))
  
  
  #  initialize matrix
  taxmat <- matrix("", nrow = length(col_names), ncol = length(prefixes),
                   dimnames = list(NULL, prefixes))
  
  # Fill in the matrix
  for (i in 1:length(col_names)) {
    for (j in 1:length(split_values[[i]])) {
      prefix <- gsub("__.*", "", split_values[[i]][j])
      value <- gsub("^[^_]+__", "", split_values[[i]][j])
      taxmat[i, prefix] <- value
    }
  }
  
  #remove extra columns
  cols_to_remove <- c("_Incertae_Sedis", "Thermus", "_Incertae_sedis")
  taxmat <- taxmat[, !(colnames(taxmat) %in% cols_to_remove)]
  
  #rename rows and columns
  new_colnames <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  colnames(taxmat) <- new_colnames
  rownames(taxmat) <- paste0("OTU", 1:num_samples)
  
  # initialize OTU matrix
  otumat <- t(kraken_data)
  rownames(otumat) <- paste0("OTU", 1:num_samples)
  
  OTU = otu_table(otumat, taxa_are_rows = TRUE) # OTU matrix
  TAX = tax_table(taxmat) # taxonomy matrix
  sampledata = sample_data(as.data.frame(kraken_meta)) #add metadata
  
  # create phyloseq object
  physeq = phyloseq(OTU, TAX, sampledata)
  
  return(physeq)
  
}