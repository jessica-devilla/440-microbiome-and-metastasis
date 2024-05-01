#Importing Data
# Clean environment -------------------------------------------------------
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# Print a starting message
cat("Importing data from Poore et al...\n")


# load the libraries
suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(dplyr)
  library(corrr)
  library(ggcorrplot)
  library(FactoMineR)
  library(devtools)
  library(ggbiplot)
  library(factoextra)
  library(ggrepel)
})


# import data -------------------------------------------------------------

##### IMPORT DATA AND FORMAT

cat("Downloading files from github...\n")
# import kraken data from poore et al
kraken_url <- "https://media.githubusercontent.com/media/jessica-devilla/JD_20_440_pset6/main/data/Kraken-TCGA-Voom-SNM-Plate-Center-Filtering-Data.csv"
kraken_data <- read_csv(url(kraken_url),show_col_types = FALSE)
kraken_df <- as.data.frame(kraken_data, stringsAsFactors = FALSE)

# import kraken metadata from poore et al
kraken_meta_url <- "https://media.githubusercontent.com/media/jessica-devilla/JD_20_440_pset6/main/data/Metadata-TCGA-Kraken-17625-Samples.csv"
kraken_metadata <-read_csv(url(kraken_meta_url),show_col_types = FALSE)
kraken_metadata_df <- as.data.frame(kraken_metadata, stringsAsFactors = FALSE)

cat("Saving RDS files...\n")
# Save the dataframes as an R file
saveRDS(kraken_df, file = "data/kraken_df.RDS")
saveRDS(kraken_df, file = "data/kraken_metadatadf.RDS")

# Subset the dataframe to get only values from COAD patients
kraken_metaCOAD <- subset(kraken_metadata_df, disease_type == "Colon Adenocarcinoma" & sample_type =='Primary Tumor' & pathologic_stage_label != "Not available")
#dim(kraken_metaCOAD)
#get the patient ids from COAD metadata
ids <- kraken_metaCOAD[,1]
kraken_COAD <- kraken_df[kraken_df[,1] %in% ids,]

#Group Stages Together
kraken_metaCOAD$pathologic_stage_label <- gsub("Stage IV([A-C])?", "Stage IV", kraken_metaCOAD$pathologic_stage_label)
kraken_metaCOAD$pathologic_stage_label <- gsub("Stage III([A-C])?", "Stage III", kraken_metaCOAD$pathologic_stage_label)
kraken_metaCOAD$pathologic_stage_label <- gsub("Stage II([A-C])?", "Stage II", kraken_metaCOAD$pathologic_stage_label)
kraken_metaCOAD$pathologic_stage_label <- gsub("Stage I([A-C])?", "Stage I", kraken_metaCOAD$pathologic_stage_label)

#Filter Taxonomy Unit into Genus Labels
kraken_COAD_genus <- kraken_COAD
colnames(kraken_COAD_genus) <- sub(".*__(.*)$", "\\1", colnames(kraken_COAD_genus))

#Filtering Metadata by Submitting Center and Experimental Type - Metadata and COAD data
split_metadata_submittingcenter <- split(kraken_metaCOAD, kraken_metaCOAD$data_submitting_center_label)
unique_centers <- unique(kraken_metaCOAD$data_submitting_center_label)
unique_source <- unique(kraken_metaCOAD$tissue_source_site_label)

for (center in unique_centers) {
  dataset_name <- paste0("kraken_meta_", gsub(" ", "", center)) # Generate a unique variable name
  assign(dataset_name, split_metadata_submittingcenter[[center]]) # Assign the dataset to the variable
  
  current_dataset <- get(dataset_name)

  data_kraken_name <- paste0("kraken_COADg_", gsub(" ", "", center)) # Generate a unique variable name for kraken_COAD
  filtered_kraken <- kraken_COAD_genus %>%
    filter(...1 %in% current_dataset$...1) # Replace ...1 with the appropriate column name from kraken_COAD
  
  assign(data_kraken_name, filtered_kraken) # Assign the filtered kraken_COAD dataset
  split_by_strategy <- split(current_dataset, current_dataset$experimental_strategy)

  for (strategy in unique(current_dataset$experimental_strategy)) {
    strategy_name <- paste0(dataset_name, "_", gsub(" ", "", strategy))
    assign(strategy_name, split_by_strategy[[strategy]])
    print(strategy_name)
    current_strategy <- get(strategy_name)
    
    data_strategy_name <- paste0("kraken_COADg_", gsub(" ", "", center), "_", gsub(" ", "", strategy))
    filtered_kraken_strategy <- kraken_COAD_genus %>%
      filter(...1 %in%current_strategy$...1) # Replace ...1 with the appropriate column name from kraken_COAD
    assign(data_strategy_name, filtered_kraken_strategy) # Assign the filtered kraken_COAD dataset for the current strategy
  }
}

#Filtering by Tissue Source Site
split_meta_source <- split(kraken_metaCOAD, kraken_metaCOAD$tissue_source_site_label)

for (source in unique_source){
  dataset_name <- paste0("kraken_meta_", gsub(" ", "", source)) # Generate a unique variable name
  assign(dataset_name, split_meta_source[[source]]) # Assign the dataset to the variable
  
  current_dataset <- get(dataset_name)
  
  data_kraken_name <- paste0("kraken_COADg_", gsub(" ", "", source)) # Generate a unique variable name for kraken_COAD
  filtered_kraken <- kraken_COAD_genus %>%
    filter(...1 %in% current_dataset$...1) # Replace ...1 with the appropriate column name from kraken_COAD
  
  assign(data_kraken_name, filtered_kraken) # Assign the filtered kraken_COAD dataset
}










#Clean Dataframe for only RNA-Seq
kraken_metaCOAD_RNASeq <- subset(kraken_metaCOAD, experimental_strategy == 'RNA-Seq')
row.names(kraken_metaCOAD_RNASeq) <- kraken_metaCOAD_RNASeq$...1
kraken_COADdata_RNASeq <- kraken_COAD %>%
  filter(...1 %in% kraken_metaCOAD_RNASeq$...1)
row.names(kraken_COADdata_RNASeq) <- kraken_COADdata_RNASeq$...1
kraken_COADRNA_clean <- subset(kraken_COADdata_RNASeq, select = -c(...1))
kraken_metaRNA_clean <- subset(kraken_metaCOAD_RNASeq, select = -c(...1))


#Filter for IlluminaGA
kraken_metaCOAD_RNASeq_IlluminaGA <- subset(kraken_metaCOAD_RNASeq, platform == "Illumina GA")
kraken_COADdata_RNASeq_IlluminaGA <- kraken_COADdata_RNASeq %>%
  filter(row.names(.) %in% row.names(kraken_metaCOAD_RNASeq_IlluminaGA))
kraken_COADRNA_Illumina_clean <- subset(kraken_COADdata_RNASeq_IlluminaGA, select = -c(...1))
kraken_metaCOADRNA_Illumina_clean <- subset(kraken_metaCOAD_RNASeq_IlluminaGA, select = -c(...1))

#Filter for UNC
kraken_metaCOAD_RNASeq_IlluminaGA_UNC <- subset(kraken_metaCOAD_RNASeq_IlluminaGA, data_submitting_center_label == "University of North Carolina")
kraken_COADdata_RNASeq_IlluminaGA_UNC <- kraken_COADdata_RNASeq_IlluminaGA %>%
  filter(row.names(.) %in% row.names(kraken_metaCOAD_RNASeq_IlluminaGA_UNC))
kraken_COADRNA_Illumina_UNC_clean <- subset(kraken_COADdata_RNASeq_IlluminaGA_UNC, select = -c(...1))
kraken_metaCOADRNA_Illumina_UNC_clean <- subset(kraken_metaCOAD_RNASeq_IlluminaGA_UNC, select = -c(...1))

#Filtering/Importing Non-Normalized Data,
kraken_orig_otu <- Kraken_TCGA_Raw_Data_17625_Samples
kraken_orig_otu_df <- as.data.frame(kraken_orig_otu)
kraken_raw_COADRNA_IlluminaGA_UNC <- kraken_orig_otu_df %>%
  filter(...1 %in% kraken_metaCOAD_RNASeq_IlluminaGA_UNC$...1)
row.names(kraken_raw_COADRNA_IlluminaGA_UNC) <- kraken_raw_COADRNA_IlluminaGA_UNC$...1
kraken_raw_COADRNA_IlluminaGA_UNC_clean <- subset(kraken_raw_COADRNA_IlluminaGA_UNC, select = -c(...1))


