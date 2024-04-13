#Different Data Filtering Techniques
kraken_orig_otu <- Kraken_TCGA_Raw_Data_17625_Samples
kraken_orig_otu_df <- as.data.frame(kraken_orig_otu)

# Subset the dataframe to get only values from COAD patients
kraken_metaCOAD <- subset(kraken_metadata_df, disease_type == "Colon Adenocarcinoma")
ids <- kraken_metaCOAD[,1]
kraken_COAD <- kraken_df[kraken_df[,1] %in% ids,]


#Filter Meta-Data by RNA-Seq Type - This is the Pre-Normalized Data from Poore et. al
kraken_metaCOAD_RNASeq <- subset(kraken_metaCOAD, sample_type == "Primary Tumor" & experimental_strategy == 'RNA-Seq' & pathologic_stage_label != "Not available")
kraken_metaCOAD_WGS <- subset(kraken_metaCOAD, sample_type == "Primary Tumor" & experimental_strategy == 'WGS' & pathologic_stage_label != "Not available")
kraken_COADdata_RNASeq <- kraken_COAD %>%
  filter(...1 %in% kraken_metaCOAD_RNASeq$...1)
row.names(kraken_COADdata_RNASeq) <- kraken_COADdata_RNASeq$...1
kraken_COADdata_RNASeq$...1 <- NULL
kraken_metaCOAD_RNASeq$pathologic_stage_label <- gsub("Stage IV([A-C])?", "Stage IV", kraken_metaCOAD_RNASeq$pathologic_stage_label)
kraken_metaCOAD_RNASeq$pathologic_stage_label <- gsub("Stage III([A-C])?", "Stage III", kraken_metaCOAD_RNASeq$pathologic_stage_label)
kraken_metaCOAD_RNASeq$pathologic_stage_label <- gsub("Stage II([A-C])?", "Stage II", kraken_metaCOAD_RNASeq$pathologic_stage_label)
kraken_metaCOAD_RNASeq$pathologic_stage_label <- gsub("Stage I([A-C])?", "Stage I", kraken_metaCOAD_RNASeq$pathologic_stage_label)
row.names(kraken_metaCOAD_RNASeq) <- kraken_metaCOAD_RNASeq$...1
kraken_metaCOAD_RNASeq_filtered <- kraken_metaCOAD_RNASeq
kraken_metaCOAD_RNASeq_filtered$...1 <- NULL

kraken_metaCOAD_RNASeq_IlluminaGA <- subset(kraken_metaCOAD_RNASeq, platform == "Illumina GA")
kraken_COADdata_RNASeq_IlluminaGA <- kraken_COADdata_RNASeq %>%
  filter(row.names(.) %in% row.names(kraken_metaCOAD_RNASeq_IlluminaGA))

#Filter Orig Data (Pre-Normalization) by RNA-Seq Types
kraken_nonorm_COADdata_RNASeq <- kraken_orig_otu_df %>%
  filter(...1 %in% kraken_metaCOAD_RNASeq$...1)
row.names(kraken_nonorm_COADdata_RNASeq) <- kraken_nonorm_COADdata_RNASeq$...1
kraken_nonorm_COADdata_RNASeq$...1 <- NULL

kraken_COADdata_nnRNASeq_IlluminaGA <- kraken_nonorm_COADdata_RNASeq %>%
  filter(row.names(.) %in% row.names(kraken_metaCOAD_RNASeq_IlluminaGA))

