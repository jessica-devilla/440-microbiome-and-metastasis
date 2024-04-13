if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install("MMUPHin")
BiocManager::install("ccrepe")
library(ccrepe)
merged_data_COAD <- merge(kraken_metaCOAD_subset, kraken_COAD, by = "...1")

subset_and_remove_column <- function(df, stage_labels) {
  subset_data <- subset(df, pathologic_stage_label %in% stage_labels)
  subset_data <- subset_data[, !names(subset_data) %in% "pathologic_stage_label"]
  subset_data <- subset_data[, !names(subset_data) %in% "...1"]
  rownames(subset_data) <- NULL
  return(subset_data)
}
stageI_COADdata <- subset_and_remove_column(merged_data_COAD, c("Stage IA", "Stage IB", "Stage I"))
stageII_COADdata <- subset_and_remove_column(merged_data_COAD, c("Stage IIA", "Stage IIB", "Stage II"))
stageIII_COADdata <- subset_and_remove_column(merged_data_COAD, c("Stage IIIA", "Stage IIIB", "Stage III"))
stageIV_COADdata <- subset_and_remove_column(merged_data_COAD, c("Stage IVA", "Stage IVB", "Stage IV"))

install.packages('infotheo')


ccrepe_function <- ccrepe(x = stageII_COADdata, y = stageIII_COADdata, sim.score = nc.score(x=stageII_COADdata, y = stageIII_COADdata), iterations = 20, min.subj = 10)
