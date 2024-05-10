#PCA FOR TISSUE SOURCE SITE
colnames(kraken_COADdata_RNASeq_IlluminaGA_UNC)[1] <- "id"
colnames(kraken_metaCOAD_RNASeq_IlluminaGA_UNC)[1] <- "id"
kraken_merge_RNAIlluminaUNC <- merge(kraken_COADdata_RNASeq_IlluminaGA_UNC,kraken_metaCOAD_RNASeq_IlluminaGA_UNC[,c('id','pathologic_stage_label', 'tissue_source_site_label')], by='id',all.x=TRUE)
contaminant_columns <- grepl("contaminant", names(kraken_merge_RNAIlluminaUNC), ignore.case = TRUE)
# Subset the dataframe to exclude contaminant columns
kraken_merge_RNAIlluminaUNC <- kraken_merge_RNAIlluminaUNC[, !contaminant_columns]
colnames(kraken_merge_RNAIlluminaUNC) <- sub(".*__(.*)$", "\\1", colnames(kraken_merge_RNAIlluminaUNC))
kraken_RNA_Ill_UNC <- kraken_merge_RNAIlluminaUNC[ , !(names(kraken_merge_RNAIlluminaUNC) %in% c("id", "pathologic_stage_label", 'tissue_source_site_label'))]
kraken_RNA_Ill_UNC_pca <- scale(kraken_RNA_Ill_UNC)
# PCA
data2.pca <- prcomp(kraken_RNA_Ill_UNC_pca, center = FALSE, scale = FALSE)
data2.pca$rotation[1:10, 1:2]
pdf("pca_UNC_RNA_illumina_filtered.pdf")  # Open a PDF graphics device and specify the filename
fviz_eig(data2.pca, addlabels = TRUE)  # Generate the plot
dev.off()
#PLOTTING PCA WITH ELIPSES 
library(factoextra)
library(ggplot2)
library(RColorBrewer)
data2.pca <- prcomp(kraken_RNA_Ill_UNC_pca, center = FALSE, scale = FALSE)
pca_cos2_plot <- fviz_pca_ind(data2.pca, col.ind="cos2") +
  scale_color_gradient2(low="white", mid="blue", high="red", midpoint=0.6) +
  stat_ellipse(aes(color = kraken_merge_RNAIlluminaUNC$tissue_source_site_label), type = "norm")
# Plot PCA colored by group with confidence ellipses
pca_group_plot <- fviz_pca_ind(data2.pca, label = "none", col.ind = kraken_merge_RNAIlluminaUNC$tissue_source_site_label) +
  stat_ellipse(aes(color = kraken_merge_RNAIlluminaUNC$tissue_source_site_label), type = "norm")
pdf("pca_krakenUNCRNA_location_ovalelipse.pdf")
print(pca_cos2_plot)
print(pca_group_plot)
dev.off()
#PLOTTING WITH SHADED ELIPSES
pca_data <- data.frame(data2.pca$x[, 1:2])  # Extract the first two principal components
#Add Grouping Variable
pca_data$group <- kraken_merge_RNAIlluminaUNC$tissue_source_site_label
color_map <- brewer.pal(n = length(unique(pca_data$group)), name = "Set1")
pca_data$outline_color <- color_map[as.numeric(factor(pca_data$group))]
pca_cos2_plot <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
  stat_ellipse(aes(group = group, color = outline_color, fill = outline_color), type = "norm", geom = "polygon", alpha = 0.2) +
  scale_color_manual(values = color_map, guide = "none") +  # Specify colors
  scale_fill_manual(values = color_map, name = "Location", labels = unique(pca_data$group) ) +  # Specify legend labels
  theme_minimal()
# Save plot to PDF
pdf("pca_krakenUNCRNA_elipse_location.pdf")
print(pca_cos2_plot)
dev.off()


#PCA FOR SUBMITTING CENTER
colnames(kraken_COAD_genus)[1] <- "id"
colnames(kraken_meta_COAD_genus)[1] <- "id"
kraken_merge <- merge(kraken_COAD_genus,kraken_meta_COAD_genus[,c('id','pathologic_stage_label', 'data_submitting_center_label')], by='id',all.x=TRUE)
contaminant_columns <- grepl("contaminant", names(kraken_merge), ignore.case = TRUE)
kraken_merge <- kraken_merge[, !contaminant_columns]
krakenpca_total <- kraken_merge[ , !(names(kraken_merge) %in% c("id", "pathologic_stage_label", 'data_submitting_center_label'))]
krakenpca_total <- scale(krakenpca_total)
#PCA
data2.pca <- prcomp(krakenpca_total, center = FALSE, scale = FALSE)
pca_data <- data.frame(data2.pca$x[, 1:2])  # Extract the first two principal components
pca_data$group <- kraken_merge$data_submitting_center_label
color_map <- brewer.pal(n = length(unique(pca_data$group)), name = "Set1")
pca_data$outline_color <- color_map[as.numeric(factor(pca_data$group))]

pca_cos2_plot <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
  stat_ellipse(aes(group = group, color = outline_color, fill = outline_color), type = "norm", geom = "polygon", alpha = 0.2) +
  scale_color_manual(values = color_map, guide = "none") +  # Specify colors
  scale_fill_manual(values = color_map, name = "Location", labels = unique(pca_data$group) ) +  # Specify legend labels
  theme_minimal()

# Save plot to PDF
pdf("pca_krakendatasubmit_elipse.pdf")
print(pca_cos2_plot)
dev.off()


