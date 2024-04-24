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
kraken_metaCOAD <- subset(kraken_metadata_df, disease_type == "Colon Adenocarcinoma")
#dim(kraken_metaCOAD)
#get the patient ids from COAD metadata
ids <- kraken_metaCOAD[,1]
kraken_COAD <- kraken_df[kraken_df[,1] %in% ids,]

# import the R files
kraken_COAD <- readRDS('data/kraken_COAD.RDS')
kraken_metaCOAD  <- readRDS('data/kraken_metaCOAD.RDS')



#Clean Dataframe for only RNA-Seq + Primary Tumor
kraken_metaCOAD_RNASeq <- subset(kraken_metaCOAD, sample_type == "Primary Tumor" & experimental_strategy == 'RNA-Seq' & pathologic_stage_label != "Not available")
row.names(kraken_metaCOAD_RNASeq) <- kraken_metaCOAD_RNASeq$...1
kraken_COADdata_RNASeq <- kraken_COAD %>%
  filter(...1 %in% kraken_metaCOAD_RNASeq$...1)
row.names(kraken_COADdata_RNASeq) <- kraken_COADdata_RNASeq$...1
kraken_COADRNA_clean <- subset(kraken_COADdata_RNASeq, select = -c(...1))
#Clean Stage Labels
kraken_metaCOAD_RNASeq$pathologic_stage_label <- gsub("Stage IV([A-C])?", "Stage IV", kraken_metaCOAD_RNASeq$pathologic_stage_label)
kraken_metaCOAD_RNASeq$pathologic_stage_label <- gsub("Stage III([A-C])?", "Stage III", kraken_metaCOAD_RNASeq$pathologic_stage_label)
kraken_metaCOAD_RNASeq$pathologic_stage_label <- gsub("Stage II([A-C])?", "Stage II", kraken_metaCOAD_RNASeq$pathologic_stage_label)
kraken_metaCOAD_RNASeq$pathologic_stage_label <- gsub("Stage I([A-C])?", "Stage I", kraken_metaCOAD_RNASeq$pathologic_stage_label)
kraken_metaRNA_clean <- subset(kraken_metaCOAD_RNASeq, select = -c(...1))
#Filter for IlluminaGA
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

saveRDS(kraken_COADRNA_Illumina_UNC_clean , file = "data/kraken_norm_filtered.RDS")
saveRDS(kraken_metaCOADRNA_Illumina_UNC_clean , file = "data/kraken_meta_norm_filtered.RDS")

#Filtering/Importing Non-Normalized Data

kraken_orig_otu <- read_csv('data/Kraken-TCGA-Raw-Data-17625-Samples.csv',show_col_types = FALSE)
kraken_orig_otu_df <- as.data.frame(kraken_orig_otu)
kraken_raw_COADRNA_IlluminaGA_UNC <- kraken_orig_otu_df %>%
  filter(...1 %in% kraken_metaCOAD_RNASeq_IlluminaGA_UNC$...1)
row.names(kraken_raw_COADRNA_IlluminaGA_UNC) <- kraken_raw_COADRNA_IlluminaGA_UNC$...1
kraken_raw_COADRNA_IlluminaGA_UNC_clean <- subset(kraken_raw_COADRNA_IlluminaGA_UNC, select = -c(...1))


saveRDS(kraken_raw_COADRNA_IlluminaGA_UNC_clean , file = "data/kraken_raw_filtered.RDS")



