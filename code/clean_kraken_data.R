suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

remove_viruses <- function(df) {
  virus_columns <- grepl("^k__Viruses", names(df), ignore.case = TRUE)
  df <- df[, !virus_columns]
  viroid_columns <- grepl("^k__Viroids", names(df), ignore.case = TRUE)
  df <- df[, !viroid_columns]
  return(df)
}

remove_contaminants <- function(df){
  # Subset the dataframe to exclude contaminant columns
  contaminant_columns <- grepl("contaminant", names(df), ignore.case = TRUE)
  df <- df[, !contaminant_columns]
  
}

clean_kraken_data <- function(kraken_data, kraken_meta){
  
  kraken_data <- remove_contaminants(kraken_data)
  kraken_data <- remove_viruses(kraken_data)
  
  if ("...1" %in% colnames(kraken_data)) {
    print("Removing column '...1'")
    row.names(kraken_data) <- kraken_data$...1
    kraken_data <- subset(kraken_data, select = -c(...1))
    
    row.names(kraken_meta) <- kraken_meta$...1
    kraken_meta <- subset(kraken_meta, select = -c(...1))
  }
  
  
  kraken_meta$pathologic_stage_label <- gsub("Stage IV([A-C])?", "Stage IV", kraken_meta$pathologic_stage_label)
  kraken_meta$pathologic_stage_label <- gsub("Stage III([A-C])?", "Stage III", kraken_meta$pathologic_stage_label)
  kraken_meta$pathologic_stage_label <- gsub("Stage II([A-C])?", "Stage II", kraken_meta$pathologic_stage_label)
  kraken_meta$pathologic_stage_label <- gsub("Stage I([A-C])?", "Stage I", kraken_meta$pathologic_stage_label)
  
  return(list(kraken_data = kraken_data, kraken_meta = kraken_meta))
}