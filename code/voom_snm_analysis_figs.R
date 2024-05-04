rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(phyloseq); packageVersion("phyloseq")
  library(ggplot2)
  library(mia)
  library(scater)
  library(umap)
  library(microbiome)
  library(plyr)
  library(vegan)
  library(reshape2)
  library(corrr)
  library(ggcorrplot)
  library(FactoMineR)
  library(devtools)
  library(ggbiplot)
  library(factoextra)
  library(ggrepel)
})

source("code/clean_kraken_data.R")
source("code/run_norm_func.R")
source("code/make_phyloseq_obj.R")
source("code/phyloseq_beta_diversity.R")
source("code/run_zicoseq.R")

kraken_voom_snm <- readRDS("data/kraken_COAD.RDS")
kraken_meta_voom <- readRDS("data/kraken_metaCOAD.RDS")
result_voom <- clean_kraken_data(kraken_voom_snm,kraken_meta_voom)
kraken_data <- result_voom$kraken_data
kraken_meta <- result_voom$kraken_meta

kraken_pca <- kraken_data

colnames(kraken_pca) <- sub(".*__(.*)$", "\\1", colnames(kraken_pca))

kraken_pca <- scale(kraken_pca)


# Run PCA
data.pca <- prcomp(kraken_pca, center = FALSE, scale = FALSE)

#summary(data.pca)

fviz_eig(data.pca, addlabels = TRUE)

## GRAPHS OF INDIVIDUALS
# PCA colored by cos2
fviz_pca_ind(data.pca, geom="point",col.ind="cos2")+
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=0.6)

# Color individuals by group
fviz_pca_ind(data.pca, label="none",
             col.ind = kraken_meta$pathologic_stage_label) +
  theme_minimal()+
  scale_shape_manual(values=c(19,19,19,19,19))+
  scale_color_manual(values = c("blue", "#8070FE", "#EAB606","#FC4703"))+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12)) 

ggsave(path = "figures", filename = "pca_allsamples_bystage.png", bg='white',width=5.5, height=5)


##GRAPHS OF VARIABLES
fviz_pca_var(data.pca, select.var = list(contrib = 10))


g <- fviz_pca_biplot(data.pca,
                     label='var', 
                     repel = TRUE, 
                     select.var = list(contrib = 10),
                     alpha.ind = 0.75,
                     #col.ind = "gray",
                     col.var = 'black',
                     labelsize = 4,
                     arrowsize = 1,
                     alpha.var=0.35,
                     col.ind = kraken_meta$pathologic_stage_label) +
  scale_color_manual(name= "Stage",values = c("blue", "#8070FE", "#EAB606","#FC4703"))+
  theme_minimal()+
  scale_shape_manual(values=c(19,19,19,19,19))+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12)) 

print(g)
ggsave(path = "figures", filename = "pca_allsamples_bystage_biplot.png", bg='white', width=5.5, height=5)
