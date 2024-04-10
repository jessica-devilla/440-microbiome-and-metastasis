#track 
kraken_metaCOAD_subset <- subset(kraken_metaCOAD, select = c("...1", "pathologic_stage_label"))
merged_data_COAD <- merge(kraken_metaCOAD_subset, kraken_COAD, by = "...1")


#cleaningdata
subset_and_remove_column <- function(df, stage_labels) {
  subset_data <- subset(df, pathologic_stage_label %in% stage_labels)
  subset_data <- subset_data[, !names(subset_data) %in% "pathologic_stage_label"]
  return(subset_data)
}
stageI_COADdata <- subset_and_remove_column(merged_data_COAD, c("Stage IA", "Stage IB", "Stage I"))
stageII_COADdata <- subset_and_remove_column(merged_data_COAD, c("Stage IIA", "Stage IIB", "Stage II"))
stageIII_COADdata <- subset_and_remove_column(merged_data_COAD, c("Stage IIIA", "Stage IIIB", "Stage III"))
stageIV_COADdata <- subset_and_remove_column(merged_data_COAD, c("Stage IVA", "Stage IVB", "Stage IV"))

