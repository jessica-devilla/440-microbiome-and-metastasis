library(GUniFrac)
library(dplyr)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")


rm(list = ls(all.names = TRUE))

remove_viruses <- function(df) {
  virus_columns <- grepl("^k__Viruses", names(df), ignore.case = TRUE)
  df <- df[, !virus_columns]
  return(df)
}

remove_contaminants <- function(df){
  # Subset the dataframe to exclude contaminant columns
  contaminant_columns <- grepl("contaminant", names(kraken_data), ignore.case = TRUE)
  kraken_data <- kraken_data[, !contaminant_columns]
  
}

kraken_meta <- readRDS("data/kraken_meta_norm_filtered.RDS")
kraken_data<- readRDS("data/kraken_norm_filtered.RDS")

kraken_data <- remove_viruses(kraken_data)
kraken_data <- remove_contaminants(kraken_data)

