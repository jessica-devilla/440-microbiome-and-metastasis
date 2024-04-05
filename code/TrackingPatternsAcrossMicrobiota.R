#track 
kraken_metaCOAD_subset <- subset(kraken_metaCOAD, select = c("...1", "pathologic_stage_label"))
merged_data_COAD <- merge(kraken_metaCOAD_subset, kraken_COAD, by = "...1")

