kraken_meta_COAD_genus <- kraken_metaCOAD
colnames(kraken_COAD_genus)[1] <- "id"
colnames(kraken_meta_COAD_genus)[1] <- "id"
kraken_genusabundance <- merge(kraken_COAD_genus, kraken_meta_COAD_genus[c('id','pathologic_stage_label')], by='id',all.x=TRUE)
