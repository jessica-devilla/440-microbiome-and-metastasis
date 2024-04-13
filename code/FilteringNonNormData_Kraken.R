#Different Data Filtering Techniques

#Filter Meta-Data by Experimental Type - This is the Pre-Normalized Data from Poore et. al
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
kraken_metaCOAD_subset$...1 <- NULL


#Filter Orig Data (Pre-Normalization) by experimental Types
kraken_nonorm_COADdata_RNASeq <- kraken_orig_otu_df %>%
  filter(...1 %in% kraken_metaCOAD_RNASeq$...1)
row.names(kraken_nonorm_COADdata_RNASeq) <- kraken_nonorm_COADdata_RNASeq$...1
kraken_nonorm_COADdata_RNASeq$...1 <- NULL
