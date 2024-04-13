#Trying Different Normalization Techniques on the Raw Kraken Output
install.packages("devtools")  # Install devtools package if you haven't already
devtools::install_github("joey711/phyloseq")


kraken_orig_otu <- Kraken_TCGA_Raw_Data_17625_Samples
kraken_orig_otu_df <- as.data.frame(kraken_orig_otu)
kraken_data_nonorm_clean <- kraken_orig_otu_df
rownames(kraken_data_nonorm_clean) <- kraken_data_nonorm_clean$...1
kraken_data_nonorm_clean$...1 <- NULL


#TSS
total_counts <- rowSums(kraken_data_nonorm_clean)
kraken_nn_tss <- t(apply(kraken_data_nonorm_clean, 1, function(x) x / sum(x) * median(total_counts)))
total_counts_rnaseq <- rowSums(kraken_nonorm_COADdata_RNASeq)
kraken_nn_tss_rnaseq <- t(apply(kraken_nonorm_COADdata_RNASeq, 1, function(x) x / sum(x) * median(total_counts)))



# Relative Abundance
kraken_nn_rabund <- t(apply(kraken_data_nonorm_clean, 1, function(x) x / sum(x)))
kraken_nn_rabund_rnaseq <- t(apply(kraken_nonorm_COADdata_RNASeq, 1, function(x) x / sum(x)))

#Rarefaction
library(phyloseq)

rarefy_phyloseq <- function(otu_df, sample.size = 1900, rngseed = 1) {
  physeq <- phyloseq(otu_table(otu_df, taxa_are_rows = FALSE))
  physeq_rarefied <- rarefy_even_depth(physeq, rngseed = rngseed, sample.size = sample.size)
  otu_table_rarefied <- as(otu_table(physeq_rarefied), "matrix")
  return(otu_table_rarefied)
}

rarefied_krakenCOAD <- rarefy_phyloseq(kraken_data_nonorm_clean)
rarefied_krakenCOAD_RNAseq <- rarefy_phyloseq(kraken_nonorm_COADdata_RNASeq)

#Quantile Normalization
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("preprocessCore")
perform_quantile_normalization <- function(data) {
  library(preprocessCore)
  data_matrix <- as.matrix(data)
  normalized_data <- normalize.quantiles(data_matrix)
  normalized_df <- as.data.frame(normalized_data)
  return(normalized_df)
}
quantile_krakenCOAD <- perform_quantile_normalization(kraken_data_nonorm_clean)
quantile_krakenCOAD_Rnaseq <- perform_quantile_normalization(kraken_nonorm_COADdata_RNASeq)

#RAIDA Normalization
#GMPR - only on RNA-Seq Normalization

#Normalization methods towards mitigating zero-inflation