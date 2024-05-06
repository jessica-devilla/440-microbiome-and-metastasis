
#!/usr/bin/env/Rscript --vanilla

#This script imports data and metadata from Poore et al and saves r dataframes as RDS files

# Clean environment -------------------------------------------------------
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# Print a starting message
cat("Importing data from Poore et al...\n")


# load the libraries
suppressPackageStartupMessages({
  library(readr)
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

kraken_raw <- read_csv("data/Kraken-TCGA-Raw-Data-17625-Samples.csv",show_col_types = FALSE)
kraken_raw_df <- as.data.frame(kraken_raw,stringsAsFactors=FALSE)

cat("Saving RDS files...\n")
# Save the dataframes as an R file
saveRDS(kraken_df, file = "data/kraken_df.RDS")
saveRDS(kraken_metadata_df, file = "data/kraken_metadatadf.RDS")
saveRDS(kraken_raw_df, file = "data/kraken_raw.RDS")

# Subset the dataframe to get only values from COAD patients
kraken_metaCOAD <- subset(kraken_metadata_df, disease_type == "Colon Adenocarcinoma")
#dim(kraken_metaCOAD)
#get the patient ids from COAD metadata
ids <- kraken_metaCOAD[,1]
kraken_COAD <- kraken_df[kraken_df[,1] %in% ids,]
#dim(kraken_COAD) # check to see if dimensions matched
kraken_COAD_raw <- kraken_raw_df[kraken_raw_df[,1] %in% ids,]

# Save the dataframes as an R file
saveRDS(kraken_COAD, file = "data/kraken_COAD.RDS")
saveRDS(kraken_metaCOAD, file = "data/kraken_metaCOAD.RDS")
saveRDS(kraken_COAD_raw, file = "data/kraken_COAD_raw.RDS")


