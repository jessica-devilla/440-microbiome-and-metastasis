#install.packages("GUniFrac")

library(GUniFrac)
library(dplyr)

####  load data
kraken_filt_raw <- readRDS("data/kraken_raw_filtered.RDS")
kraken_filt_norm <- readRDS("data/kraken_norm_filtered.RDS")

data(throat.otu.tab)
data(throat.meta)
comm <- t(throat.otu.tab)
meta.dat <- throat.meta
meta.dat


kraken_filt$pathologic_stage_label <- gsub("Stage IV([A-C])?", "Stage IV", kraken_filt$pathologic_stage_label)
kraken_filt$pathologic_stage_label <- gsub("Stage III([A-C])?", "Stage III", kraken_filt$pathologic_stage_label)
kraken_filt$pathologic_stage_label <- gsub("Stage II([A-C])?", "Stage II", kraken_filt$pathologic_stage_label)
kraken_filt$pathologic_stage_label <- gsub("Stage I([A-C])?", "Stage I", kraken_filt$pathologic_stage_label)
row.names(kraken_filt) <- kraken_filt$...1
kraken_filt$...1 <- NULL



zicoseq_data <- t(kraken_COAD_subset)
zico_test_1 <- ZicoSeq(meta.dat = kraken_filt,  feature.dat = zicoseq_data, grp.name = 'pathologic_stage_label', adj.name = NULL, feature.dat.type = 'other', prev.filter = 0, mean.abund.filter = 0,  max.abund.filter = 0.002, min.prop = 0, 
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
