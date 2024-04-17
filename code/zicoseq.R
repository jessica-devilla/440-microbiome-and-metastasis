install.packages("GUniFrac")
library(GUniFrac)
library(dplyr)


#DataCleaning
kraken_metaCOAD_subset <- subset(kraken_metaCOAD, sample_type == "Primary Tumor" & pathologic_stage_label != "Not available")



kraken_metaCOAD_subset <- kraken_metaCOAD_subset[, c("...1", "pathologic_stage_label")]
kraken_COAD_subset <- kraken_COAD %>%
  filter(...1 %in% kraken_metaCOAD_subset$...1)

row.names(kraken_COAD_subset) <- kraken_COAD_subset$...1
kraken_COAD_subset$...1 <- NULL



kraken_metaCOAD_subset$pathologic_stage_label <- gsub("Stage IV([A-C])?", "Stage IV", kraken_metaCOAD_subset$pathologic_stage_label)
kraken_metaCOAD_subset$pathologic_stage_label <- gsub("Stage III([A-C])?", "Stage III", kraken_metaCOAD_subset$pathologic_stage_label)
kraken_metaCOAD_subset$pathologic_stage_label <- gsub("Stage II([A-C])?", "Stage II", kraken_metaCOAD_subset$pathologic_stage_label)
kraken_metaCOAD_subset$pathologic_stage_label <- gsub("Stage I([A-C])?", "Stage I", kraken_metaCOAD_subset$pathologic_stage_label)
row.names(kraken_metaCOAD_subset) <- kraken_metaCOAD_subset$...1
kraken_metaCOAD_subset$...1 <- NULL

zicoseq_data <- t(kraken_COAD_subset)
zico_test_1 <- ZicoSeq(meta.dat = kraken_metaCOAD_subset,  feature.dat = zicoseq_data, grp.name = 'pathologic_stage_label', adj.name = NULL, feature.dat.type = 'other', prev.filter = 0, mean.abund.filter = 0,  max.abund.filter = 0.002, min.prop = 0, 
                       # Winsorization to replace outliers
                       is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                       # Posterior sampling to impute zeros
                       is.post.sample = TRUE, post.sample.no = 25, 
                       # Use the square-root transformation
                       link.func = list(function (x) x^0.5), stats.combine.func = max,
                       # Permutation-based multiple testing correction
                       perm.no = 99,  strata = NULL, 
                       # Reference-based multiple stage normalization
                       ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                       # Family-wise error rate control
                       is.fwer = FALSE,
                       verbose = TRUE, return.feature.dat = TRUE)
