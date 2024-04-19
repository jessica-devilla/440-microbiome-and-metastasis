
colnames(kraken_COADdata_RNASeq_IlluminaGA_UNC)[1] <- "id"
colnames(kraken_metaCOAD_RNASeq_IlluminaGA_UNC)[1] <- "id"

kraken_merge_RNAIlluminaUNC <- merge(kraken_COADdata_RNASeq_IlluminaGA_UNC,kraken_metaCOAD_RNASeq_IlluminaGA_UNC[,c('id','pathologic_stage_label', 'tissue_source_site_label')], by='id',all.x=TRUE)
contaminant_columns <- grepl("contaminant", names(kraken_merge_RNAIlluminaUNC), ignore.case = TRUE)

# Subset the dataframe to exclude contaminant columns
kraken_merge_RNAIlluminaUNC <- kraken_merge_RNAIlluminaUNC[, !contaminant_columns]


contaminant_columns <- grepl("contaminant", names(kraken_merge_RNAIlluminaUNC), ignore.case = TRUE)


colnames(kraken_merge_RNAIlluminaUNC) <- sub(".*__(.*)$", "\\1", colnames(kraken_merge_RNAIlluminaUNC))
kraken_RNA_Ill_UNC <- kraken_merge_RNAIlluminaUNC[ , !(names(kraken_merge_RNAIlluminaUNC) %in% c("id", "pathologic_stage_label", 'tissue_source_site_label'))]
kraken_RNA_Ill_UNC_pca <- scale(kraken_RNA_Ill_UNC)

# METHOD 1
#corr_matrix <- cor(kraken_pca)
#data.pca <- princomp(corr_matrix)
#data.pca$loadings[1:10, 1:2]

# METHOD 2
data2.pca <- prcomp(kraken_RNA_Ill_UNC_pca, center = FALSE, scale = FALSE)
data2.pca$rotation[1:10, 1:2]

summary(data2.pca)
pdf("pca_UNC_RNA_illumina_filtered.pdf")  # Open a PDF graphics device and specify the filename
fviz_eig(data2.pca, addlabels = TRUE)  # Generate the plot
dev.off()


#Plotting by Location with Elipses 
library(factoextra)
library(ggplot2)
pdf("pca_krakenUNCRNA_location_ovalelipse.pdf")
# Perform PCA
data2.pca <- prcomp(kraken_RNA_Ill_UNC_pca, center = FALSE, scale = FALSE)

pca_cos2_plot <- fviz_pca_ind(data2.pca, col.ind="cos2") +
  scale_color_gradient2(low="white", mid="blue", high="red", midpoint=0.6) +
  stat_ellipse(aes(color = kraken_merge_RNAIlluminaUNC$tissue_source_site_label), type = "norm")

# Plot PCA colored by group with confidence ellipses
pca_group_plot <- fviz_pca_ind(data2.pca, label = "none", col.ind = kraken_merge_RNAIlluminaUNC$tissue_source_site_label) +
  stat_ellipse(aes(color = kraken_merge_RNAIlluminaUNC$tissue_source_site_label), type = "norm")

print(pca_cos2_plot)
print(pca_group_plot)
dev.off()



#Plotting with Shaded Elipses
library(ggplot2)
library(factoextra)
library(RColorBrewer)

# Perform PCA
data2.pca <- prcomp(kraken_RNA_Ill_UNC_pca, center = FALSE, scale = FALSE)

# Extract PCA results
pca_data <- data.frame(data2.pca$x[, 1:2])  # Extract the first two principal components

# Add grouping variable
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










#Plotting by Location with Dots
pdf("pca_krakenUNCRNA_location.pdf")

# PCA plot colored by cos2
fviz_pca_ind(data2.pca, geom="point", col.ind="cos2") +
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=0.6)

# PCA plot colored by group
fviz_pca_ind(data2.pca, label = "none", col.ind = kraken_merge_RNAIlluminaUNC$tissue_source_site_label)


# Biplot of PCA variables
pca_plot <- fviz_pca_biplot(data2.pca,
                            label='var', 
                            repel = TRUE, 
                            select.var = list(contrib = 5),
                            alpha.ind = 0.5,
                            col.ind = "gray",
                            col.var = 'blue',
                            labelsize = 3,
                            arrowsize = 1,
                            alpha.var=0.15)

print(pca_plot)
dev.off()


