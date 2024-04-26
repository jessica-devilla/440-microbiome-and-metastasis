library(GUniFrac)
library(dplyr)

#source("code/pca_by_phylum_family_genus.R")

rm(list = ls(all.names = TRUE))

remove_viruses_contams <- function(df) {
  contaminant_columns <- grepl("^k__Viruses", names(df), ignore.case = TRUE)
  df <- df[, !contaminant_columns]
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