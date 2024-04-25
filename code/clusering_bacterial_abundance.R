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
  
  
})

####  load data
kraken_merge <- readRDS('data/kraken_merge.RDS')

kraken_merge <- kraken_merge %>% arrange(stage_category)

metadata <- kraken_merge[, c("id", "stage_category")]

kraken_clust <- kraken_merge[ , !(names(kraken_merge) %in% c("id","stage_category",'sample_type'))]
kraken_clust <- as.matrix(kraken_clust)
kraken_clust <- scale(kraken_clust)


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
plot(color_branches(dend_obj, k = num_clusters))

cluster_labels <- data.frame(cluster = factor(cluster_assignments), stage = factor(metadata$stage_category))
cluster_summary <- cluster_labels %>% group_by(cluster, stage) %>% summarise(count = n())

print(cluster_summary)

# Create a bar plot
p <- ggplot(cluster_summary, aes(x = cluster, y = count, fill = stage)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Cluster Summary", x = "Cluster", y = "Count") +
  theme_minimal()
ggsave(path = "figures", filename = "hierarchical_clust_summary.png", bg='white')

pheatmap(kraken_clust, cluster_rows = dend_obj, cluster_cols = FALSE,
         show_rownames = FALSE, show_colnames = FALSE)
