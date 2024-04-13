#Trying Different Normalization Techniques on the Raw Kraken Output
install.packages("devtools")  # Install devtools package if you haven't already
devtools::install_github("joey711/phyloseq")


kraken_orig_otu <- Kraken_TCGA_Raw_Data_17625_Samples
kraken_orig_otu_df <- as.data.frame(kraken_orig_otu)


#Filter Orig Data (Pre-Normalization) by experimental Types
kraken_nonorm_COADdata_RNASeq <- kraken_orig_otu_df %>%
  filter(...1 %in% kraken_metaCOAD_RNASeq$...1)
row.names(kraken_nonorm_COADdata_RNASeq) <- kraken_nonorm_COADdata_RNASeq$...1
kraken_nonorm_COADdata_RNASeq$...1 <- NULL


rownames(kraken_orig_otu_df) <- kraken_orig_otu_df$...1
kraken_orig_otu_df$...1 <- NULL


#TSS
total_counts <- rowSums(kraken_orig_otu_df)
kraken_nn_tss <- t(apply(kraken_orig_otu_df, 1, function(x) x / sum(x) * median(total_counts)))

total_counts_rnaseq <- rowSums(kraken_nonorm_COADdata_RNASeq)
kraken_nn_tss_rnaseq <- t(apply(kraken_nonorm_COADdata_RNASeq, 1, function(x) x / sum(x) * median(total_counts)))






# Relative Abundance
otu_table_rel_abun <- t(apply(otu_df, 1, function(x) x / sum(x)))


#Rarefaction
library(phyloseq)
physeq <- phyloseq(otu_table(otu_df, taxa_are_rows = FALSE))
physeq_rarefied <- rarefy_even_depth(physeq, rngseed = 1, sample.size = 1000)
otu_table_rarefied <- as(otu_table(physeq_rarefied), "matrix")

#Centered Log-Ratio Transformation
install.packages("compositions")
library(compositions)
clr_transformed <- clr(otu_df)
head(clr_transformed)

#Total Count Normalization
otu_table_tcn <- t(apply(otu_df, 1, function(x) x / sum(x) * mean(rowSums(otu_df))))


#Z-Score Transformation
otu_table_zscore <- t(apply(otu_df, 1, function(x) (x - mean(x)) / sd(x)))



