#!/usr/bin/env/Rscript --vanilla


# Clean environment -------------------------------------------------------
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage


# Print a starting message
cat("Clustering bacterial abundance data...\n")

#BiocManager::install("ComplexHeatmap")


# load the libraries
suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(here)
  library(ComplexHeatmap)
  library(dplyr)
  library(pheatmap)
  library(dendextend)
  library(cluster)
  library(mclust)
  library(aricode)
})

#import data
source("code/clean_kraken_data.R")
kraken_voom_snm <- readRDS("data/kraken_COAD.RDS")
kraken_meta_voom <- readRDS("data/kraken_metaCOAD.RDS")
result_voom <- clean_kraken_data(kraken_voom_snm,kraken_meta_voom)
kraken_data <- result_voom$kraken_data
kraken_meta <- result_voom$kraken_meta
kraken_meta['id'] <- rownames(kraken_meta)
####  load data
kraken_merge <- readRDS('data/kraken_merge.RDS')

kraken_merge <- kraken_merge %>% arrange(stage_category)

metadata <- kraken_merge[, c("id", "stage_category")]
metadata<- merge(metadata, kraken_meta[, c("id", "data_submitting_center_label","experimental_strategy","platform","vital_status","gender")], by = "id", all.x = TRUE)

kraken_clust <- kraken_merge[ , !(names(kraken_merge) %in% c("id","stage_category",'sample_type'))]
kraken_clust <- as.matrix(kraken_clust)
kraken_clust <- scale(kraken_clust)

# draw heatmap of genus clustering
heatmap <- Heatmap(kraken_clust,
        name = "Abundance",
        cluster_rows = TRUE ,cluster_columns = FALSE,
        row_dend_width = unit(0.2, "npc"),
        column_dend_height = unit(0.2, "npc"),
        show_row_names = FALSE, show_column_names = FALSE,
        row_title="ids", column_title = "Bacterial Genuses"
)
draw(heatmap)
save(heatmap, file = "figures/genus_heatmap.png")
#ggsave(path = "figures", filename = "genus_heatmap.png", bg='white')


# identify outliers
# Convert the matrix to a vector
kraken_clust_vector <- as.vector(kraken_clust)

# Plot histogram
hist(kraken_clust_vector, main = "Distribution of kraken_clust", xlab = "Values", ylab = "Frequency")

count_greater_than_5 <- sum(kraken_clust_vector > 7)

# Output the result
print(count_greater_than_5)


# do hierarchical clustering
distance_matrix <- dist(kraken_clust, method = "euclidean")
hclust_res <- hclust(distance_matrix, method = "complete")


dend <- as.dendrogram(hclust_res)
dend_obj <- as.hclust(dend)
plot(dend)


num_clusters <- 4
cluster_assignments <- cutree(hclust_res, k = num_clusters)
plot(color_branches(dend_obj, k = num_clusters, col=c("#1E88E5", "#D81B60", "#FFC107", "#004D40")))

# group ids of different stages with clustering label
cluster_labels <- data.frame(cluster = factor(cluster_assignments), stage = factor(metadata$stage_category))
cluster_summary <- cluster_labels %>%
  group_by(cluster, stage) %>%
  dplyr::summarise(count = n())

# Print cluster summary
print(cluster_summary)

#plot cluster summary
p <- ggplot(cluster_summary, aes(x = cluster, y = count, fill = stage)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Cluster Summary", x = "Cluster", y = "Count", fill = "Stage") +
  scale_fill_manual(values = c("gray", "blue", "#8070FE", "#EAB606", "#FC4703")) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

print(p)
ggsave(path = "figures", filename = "hierarchical_clust_allsamples_summary.png", bg='white', width=5.5,height=5)


pheatmap(kraken_clust, cluster_rows = dend_obj, cluster_cols = FALSE,
         show_rownames = FALSE, show_colnames = FALSE)

## Evaluate the clustering metric

# Calculate silhouette score
silhouette_score <- silhouette(cluster_assignments, dist(kraken_clust))

# Get the mean silhouette score
mean_silhouette <- mean(silhouette_score[, "sil_width"])

print(paste("Mean Silhouette Score:", mean_silhouette))

# Load true class labels (assuming you have them)
true_labels <- as.factor(metadata$stage_category)

# Calculate ARI
ari <- adjustedRandIndex(cluster_assignments, true_labels)

print(paste("Adjusted Rand Index:", ari))

# Calculate AMI
#ami <- adjustedMutualInformation(cluster_assignments, true_labels)
#print(paste("Adjusted Mutual Information:", ami))


### see if clustering happens by other variables
num_clusters <- 4
cluster_assignments <- cutree(hclust_res, k = num_clusters)
plot(color_branches(dend_obj, k = num_clusters))

cluster_labels <- data.frame(cluster = factor(cluster_assignments), stage = factor(metadata$data_submitting_center_label))
cluster_summary <- cluster_labels %>%
  group_by(cluster, stage) %>%
  dplyr::summarise(count = n())

# Print cluster summary
print(cluster_summary)


p <- ggplot(cluster_summary, aes(x = cluster, y = count, fill = stage)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Cluster Summary", x = "Cluster", y = "Count", fill = "Submitting \n Center") +
  scale_fill_manual(values = c("gray", "blue", "#8070FE", "#EAB606", "#FC4703")) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

print(p)

### see if clustering happens by other variables
num_clusters <- 2
cluster_assignments <- cutree(hclust_res, k = num_clusters)
plot(color_branches(dend_obj, k = num_clusters))

cluster_labels <- data.frame(cluster = factor(cluster_assignments), stage = factor(metadata$gender))
cluster_summary <- cluster_labels %>%
  group_by(cluster, stage) %>%
  dplyr::summarise(count = n())

# Print cluster summary
print(cluster_summary)


p <- ggplot(cluster_summary, aes(x = cluster, y = count, fill = stage)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Cluster Summary", x = "Cluster", y = "Count", fill = "Gender") +
  scale_fill_manual(values = c("gray", "blue", "#8070FE", "#EAB606", "#FC4703")) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

print(p)
