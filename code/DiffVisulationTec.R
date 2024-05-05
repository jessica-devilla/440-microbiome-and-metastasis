#Different Data Visulation
install.packages("umap")
library(umap)
# Assuming tot_TSS, tot_MED, tot_GMPR, etc. are your normalized datasets
umap_TSS <- umap(tot_TSS)
umap_MED <- umap(tot_MED)
umap_GMPR <- umap(tot_GMPR)
umap_UQ <- umap(tot_UQ)
umap_CSS <- umap(tot_CSS)
umap_DEQ <- umap(tot_DEQ)
umap_RLE <- umap(tot_RLE)
umap_RLEpos <- umap(tot_RLEpos)
umap_logcpm <- umap(tot_logcpm)
umap_rarefy <- umap(tot_rarefy)
umap_CLR <- umap(tot_CLR)
umap_CLRpos <- umap(tot_CLRpos)
umap_df <- rbind(
  data.frame(UMAP1 = umap_TSS$layout[,1], UMAP2 = umap_TSS$layout[,2], Method = "TSS"),
  data.frame(UMAP1 = umap_MED$layout[,1], UMAP2 = umap_MED$layout[,2], Method = "MED"),
  data.frame(UMAP1 = umap_GMPR$layout[,1], UMAP2 = umap_GMPR$layout[,2], Method = "GMPR"),
  data.frame(UMAP1 = umap_UQ$layout[,1], UMAP2 = umap_UQ$layout[,2], Method = "UQ"),
  data.frame(UMAP1 = umap_CSS$layout[,1], UMAP2 = umap_CSS$layout[,2], Method = "CSS"),
  data.frame(UMAP1 = umap_DEQ$layout[,1], UMAP2 = umap_DEQ$layout[,2], Method = "DEQ"),
  data.frame(UMAP1 = umap_RLE$layout[,1], UMAP2 = umap_RLE$layout[,2], Method = "RLE"),
  data.frame(UMAP1 = umap_RLEpos$layout[,1], UMAP2 = umap_RLEpos$layout[,2], Method = "RLEpos"),
  data.frame(UMAP1 = umap_logcpm$layout[,1], UMAP2 = umap_logcpm$layout[,2], Method = "logcpm"),
  data.frame(UMAP1 = umap_rarefy$layout[,1], UMAP2 = umap_rarefy$layout[,2], Method = "rarefy"),
  data.frame(UMAP1 = umap_CLR$layout[,1], UMAP2 = umap_CLR$layout[,2], Method = "CLR"),
  data.frame(UMAP1 = umap_CLRpos$layout[,1], UMAP2 = umap_CLRpos$layout[,2], Method = "CLRpos")
)
library(ggplot2)
pdf("umap_plot_total_f.pdf")
umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Method)) +
  geom_point(size = .5) +
  theme_minimal() +
  labs(title = "UMAP Visualization of Different Normalization Methods")
print(umap_plot)
dev.off()

library(Rtsne)
library(ggplot2)

# Perform t-SNE for each normalized dataset
tsne_TSS <- Rtsne(tot_TSS, perplexity = 30)  # Adjust perplexity as needed
tsne_MED <- Rtsne(tot_MED, perplexity = 30)
tsne_GMPR <- Rtsne(tot_GMPR, perplexity = 30)
tsne_UQ <- Rtsne(tot_UQ, perplexity = 30)
tsne_CSS <- Rtsne(tot_CSS, perplexity = 30)
tsne_DEQ <- Rtsne(tot_DEQ, perplexity = 30)
tsne_RLE <- Rtsne(tot_RLE, perplexity = 30)
tsne_RLEpos <- Rtsne(tot_RLEpos, perplexity = 30)
tsne_TMM <- Rtsne(tot_TMM, perplexity = 30)
tsne_logcpm <- Rtsne(tot_logcpm, perplexity = 30)
tsne_rarefy <- Rtsne(tot_rarefy, perplexity = 30)
tsne_CLR <- Rtsne(tot_CLR, perplexity = 30)
tsne_CLRpos <- Rtsne(tot_CLRpos, perplexity = 30)


# Combine t-SNE results into a single dataframe
tsne_df <- rbind(
  data.frame(tsne_TSS$Y, Method = "TSS"),
  data.frame(tsne_MED$Y, Method = "MED"),
  data.frame(tsne_GMPR$Y, Method = "GMPR"),
  data.frame(tsne_UQ$Y, Method = "UQ"),
  data.frame(tsne_CSS$Y, Method = "CSS"),
  data.frame(tsne_DEQ$Y, Method = "DEQ"),
  data.frame(tsne_RLE$Y, Method = "RLE"),
  data.frame(tsne_RLEpos$Y, Method = "RLEpos"),
  data.frame(tsne_TMM$Y, Method = "TMM"),
  data.frame(tsne_logcpm$Y, Method = "logcpm"),
  data.frame(tsne_rarefy$Y, Method = "RARE"),
  data.frame(tsne_CLR$Y, Method = "CLR"),
  data.frame(tsne_CLRpos$Y, Method = "CLRpos")
)

# Plot t-SNE visualization
tsne_plot <- ggplot(tsne_df, aes(x = X1, y = X2, color = Method)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  labs(title = "t-SNE Visualization of Different Normalization Methods")

# Save plot as PDF
pdf("tsne_plot_total_f.pdf")
print(tsne_plot)
dev.off()



# Apply PCA to each normalized dataset
pca_TSS <- prcomp(tot_TSS)
pca_MED <- prcomp(tot_MED)
pca_GMPR <- prcomp(tot_GMPR)
pca_UQ <- prcomp(tot_UQ)
pca_CSS <- prcomp(tot_CSS)
pca_DEQ <- prcomp(tot_DEQ)
pca_RLE <- prcomp(tot_RLE)
pca_RLEpos <- prcomp(tot_RLEpos)
pca_logcpm <- prcomp(tot_logcpm)
pca_rarefy <- prcomp(tot_rarefy)
pca_CLR <- prcomp(tot_CLR)
pca_CLRpos <- prcomp(tot_CLRpos)

# Combine principal components and method labels into a single dataframe
pca_df <- rbind(
  data.frame(PC1 = pca_TSS$x[,1], PC2 = pca_TSS$x[,2], Method = "TSS"),
  data.frame(PC1 = pca_MED$x[,1], PC2 = pca_MED$x[,2], Method = "MED"),
  data.frame(PC1 = pca_GMPR$x[,1], PC2 = pca_GMPR$x[,2], Method = "GMPR"),
  data.frame(PC1 = pca_UQ$x[,1], PC2 = pca_UQ$x[,2], Method = "UQ"),
  data.frame(PC1 = pca_CSS$x[,1], PC2 = pca_CSS$x[,2], Method = "CSS"),
  data.frame(PC1 = pca_DEQ$x[,1], PC2 = pca_DEQ$x[,2], Method = "DEQ"),
  data.frame(PC1 = pca_RLE$x[,1], PC2 = pca_RLE$x[,2], Method = "RLE"),
  data.frame(PC1 = pca_RLEpos$x[,1], PC2 = pca_RLEpos$x[,2], Method = "RLEpos"),
  data.frame(PC1 = pca_logcpm$x[,1], PC2 = pca_logcpm$x[,2], Method = "logcpm"),
  data.frame(PC1 = pca_rarefy$x[,1], PC2 = pca_rarefy$x[,2], Method = "rarefy"),
  data.frame(PC1 = pca_CLR$x[,1], PC2 = pca_CLR$x[,2], Method = "CLR"),
  data.frame(PC1 = pca_CLRpos$x[,1], PC2 = pca_CLRpos$x[,2], Method = "CLRpos")
)

# Visualize the PCA results colored by the normalization method
library(ggplot2)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Method)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA Visualization of Different Normalization Methods",
       color = "Normalization Method")
pdf("pca_plot_total_f.pdf")
print(pca_plot)
dev.off()




