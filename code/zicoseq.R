#install.packages("GUniFrac")

library(GUniFrac)
library(dplyr)
library(ggplot2)

#source("code/pca_by_phylum_family_genus.R")

rm(list = ls(all.names = TRUE))

remove_viruses <- function(df) {
  virus_columns <- grepl("^k__Viruses", names(df), ignore.case = TRUE)
  df <- df[, !virus_columns]
  return(df)
}

remove_contaminants <- function(df){
  # Subset the dataframe to exclude contaminant columns
  contaminant_columns <- grepl("contaminant", names(df), ignore.case = TRUE)
  df <-df[, !contaminant_columns]
  
}


run_zicoseq <- function(kraken_data, kraken_meta){
  
  # exclude viruses
  kraken_data <- remove_viruses(kraken_data)
  kraken_data <- remove_contaminants(kraken_data)
  
  #subset dataframes for comparison
  kraken_meta_sub <- kraken_meta[kraken_meta$pathologic_stage_label %in% c("Stage I", "Stage IV"), ]
  kraken_subset <- kraken_data[rownames(kraken_data) %in% rownames(kraken_meta_sub), ]
  
  
  # examine dataset
  kraken_mat <- as.matrix(kraken_subset)
  kraken_mat <- as.numeric(kraken_mat)
  print(typeof(kraken_mat))
  zero_count <- sum(kraken_subset == 0) 
  print(kraken_subset[kraken_subset==0])
  nan_count <- sum(is.nan(kraken_mat))
  print(kraken_subset[is.nan(kraken_mat)])
  
  
  #remove columns containing all zeroes from dataframe
  kraken_subset <- kraken_subset[, colSums(kraken_subset != 0) > 0]
  
  #remove columns containing less than 3 true values
  non_zero_counts <- colSums(kraken_subset != 0)
  kraken_subset <- kraken_subset[, non_zero_counts >= 3]
  
  #get only genus labels
  colnames(kraken_subset) <- sub(".*__(.*)$", "\\1", colnames(kraken_subset))
  
  zicoseq_data <- t(kraken_subset)
  #row_names <- 1:nrow(zicoseq_data)
  #rownames(zicoseq_data) <- row_names
  
  zicoObj <- ZicoSeq(meta.dat = kraken_meta_sub,  feature.dat = zicoseq_data, grp.name = 'pathologic_stage_label', 
                     adj.name = NULL, feature.dat.type = "other", prev.filter = 0, mean.abund.filter = 0,  max.abund.filter = 0, min.prop = 0, 
                     # Winsorization to replace outliers
                     is.winsor = TRUE, outlier.pct = 0.1, winsor.end = 'top',
                     # Posterior sampling 
                     is.post.sample = TRUE, post.sample.no = 25, 
                     # Use the identity function
                     link.func = list(function (x) x), stats.combine.func = max,
                     # Permutation-based multiple testing correction
                     perm.no = 99,  strata = NULL, 
                     # Reference-based multiple stage normalization
                     ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                     # Family-wise error rate control
                     is.fwer = TRUE,
                     verbose = TRUE, return.feature.dat = TRUE)
  
  
  return(list(zicoObj = zicoObj, zicoseq_data = zicoseq_data))
  
}


### RUN ZICOSEQ ON UNC DATA

kraken_meta_UNC <- readRDS("data/kraken_meta_norm_filtered.RDS")
kraken_data_UNC<- readRDS("data/kraken_norm_filtered.RDS")
result <- run_zicoseq(kraken_data_UNC, kraken_meta_UNC)

zicoObj <- result$zicoObj
zicoseq_data <- result$zicoseq_data

zico_plot_UNC <- ZicoSeq.plot(zicoObj, pvalue.type = 'p.adj.fdr', cutoff = 0.1, text.size = 10,
             out.dir = NULL, width = 10, height = 6)
print(zico_plot_UNC)

ggsave("figures/zico_plot_UNC_stagei_vs_stageiv.png", plot = zico_plot_UNC, width = 10, height = 6)

### RUN ZICOSEQ ON ALL NORM DATA

kraken_meta <- readRDS("data/kraken_metaCOAD.RDS")
kraken_data <- readRDS("data/kraken_COAD.RDS")

row.names(kraken_data) <- kraken_data$...1
kraken_data <- subset(kraken_data, select = -c(...1))

row.names(kraken_meta) <- kraken_meta$...1
kraken_meta <- subset(kraken_meta, select = -c(...1))

kraken_meta$pathologic_stage_label <- gsub("Stage IV([A-C])?", "Stage IV", kraken_meta$pathologic_stage_label)
kraken_meta$pathologic_stage_label <- gsub("Stage III([A-C])?", "Stage III", kraken_meta$pathologic_stage_label)
kraken_meta$pathologic_stage_label <- gsub("Stage II([A-C])?", "Stage II", kraken_meta$pathologic_stage_label)
kraken_meta$pathologic_stage_label <- gsub("Stage I([A-C])?", "Stage I", kraken_meta$pathologic_stage_label)

result <- run_zicoseq(kraken_data, kraken_meta,dat_type ="other", outlier_pct = 0.1)

zicoObj <- result$zicoObj
zicoseq_data <- result$zicoseq_data

zico_plot <- ZicoSeq.plot(zicoObj, pvalue.type = 'p.adj.fdr', cutoff = 0.1, text.size = 10,
                              out.dir = NULL, width = 10, height = 6)
print(zico_plot)

ggsave("figures/zico_plot_alldata_stagei_vs_stageiv.png", plot = zico_plot, width = 10, height = 6)

## RUN ZICOSEQ ON RAW DATA - STILL DEBUGGING!!

source("code/clean_kraken_data.R")


# import and process datasets of interest
kraken_meta <- readRDS("data/kraken_metaCOAD.RDS")
kraken_data <- readRDS("data/kraken_COAD_raw.RDS")

result <- clean_kraken_data(kraken_data, kraken_meta)
kraken_data <- result$kraken_data
kraken_meta <- result$kraken_meta

result <- run_zicoseq(kraken_data, kraken_meta)

zicoObj <- result$zicoObj
zicoseq_data <- result$zicoseq_data

zico_plot <- ZicoSeq.plot(zicoObj, pvalue.type = 'p.adj.fdr', cutoff = 0.1, text.size = 10,
                          out.dir = NULL, width = 10, height = 6)
print(zico_plot)

ggsave("figures/zico_plot_alldata_stagei_vs_stageiv.png", plot = zico_plot, width = 10, height = 6)
