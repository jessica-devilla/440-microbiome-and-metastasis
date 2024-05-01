kraken_meta_COAD_genus <- kraken_metaCOAD
colnames(kraken_COAD_genus)[1] <- "id"
colnames(kraken_meta_COAD_genus)[1] <- "id"
kraken_genusabundance <- merge(kraken_COAD_genus, kraken_meta_COAD_genus[c('id','pathologic_stage_label')], by='id',all.x=TRUE)
kraken_genusabundance <- kraken_genusabundance[, c("id", "pathologic_stage_label", setdiff(names(kraken_genusabundance), c("id", "pathologic_stage_label")))]
kraken_genusabundance <- kraken_genusabundance[, -1]


#First Analysis: Kraken_Total

kruskal_results <- lapply(kraken_genusabundance[, -1], function(x) {
  kruskal_result <- kruskal.test(x ~ pathologic_stage_label, data = kraken_genusabundance)
  return(kruskal_result$p.value)
})

kraken_genusabundance$pathologic_stage_label <- as.numeric(as.factor(kraken_genusabundance$pathologic_stage_label))

# Step 4: Identify taxa with significant differences between stages (low p-value) or high absolute Spearman correlation
kruskal_results_vector <- unlist(kruskal_results)
min_p_value_index <- which.min(kruskal_results_vector)
significant_taxon <- names(kruskal_results_vector[min_p_value_index])
print(significant_taxon)

strong_correlation_taxon <- names(which.max(abs(spearman_correlation)))
print(strong_correlation_taxon)


#Top 5 Significant Differences between stages and Spearman Correlation
kruskal_results_df <- data.frame(taxon = names(kruskal_results), p_value = unlist(kruskal_results))
kruskal_minp_top5 <- head(kruskal_results_df[order(kruskal_results_df$p_value), ], 5)
print(kruskal_minp_top5)
absolute_correlation <- abs(spearman_correlation)
top_5_indices_spear <- order(absolute_correlation, decreasing = TRUE)[1:5]
spear_maxabs_top5_taxa <- names(spearman_correlation)[top_5_indices_spear]
spear_maxabs_top5 <- absolute_correlation[top_5_indices_spear]

min10_kruskal_p <- sorted_p_values_kruskal[1:10]
min10_kruskal <- kruskal_results_df$taxon[order(kruskal_results_df$p_value)][1:10]
pdf("min10_pvalues_kruskal.pdf", width = 12, height = 8)
par(mar = c(10, 5, 4, 2))
plot(1:10, min10_kruskal_p, type = "b", xaxt = "n", xlab = "Taxon", ylab = "p-value", main = "Top 10 Minimum Kruskal-Wallis Test p-values")
axis(1, at = 1:10, labels = min10_kruskal, las = 2, cex.axis = 0.8)
text(1:10, par("usr")[3] - 0.05, labels = min10_kruskal, srt = 45, adj = 1, xpd = TRUE, cex = 0.8)  # Diagonal labels
dev.off()

library(dplyr)
library(ggplot2)

#Min 5 P-Value Kruskal Abundance Plot
min5_kruskal_p <- sorted_p_values_kruskal[1:5]
min5_kruskal_taxa <- kruskal_results_df$taxon[order(kruskal_results_df$p_value)][1:5]

library(dplyr)
library(ggplot2)
library(tidyr)

filter_and_plot_abundance <- function(input_df, min5_kruskal_taxa, output_file) {
  # Filter the input dataframe to include only the genus' in min5_kruskal_taxa
  filtered_df <- input_df %>%
    select(pathologic_stage_label, all_of(min5_kruskal_taxa))
  
  # Reshape the dataframe to long format for easy plotting
  filtered_df_long <- filtered_df %>%
    pivot_longer(cols = -pathologic_stage_label, names_to = "Genus", values_to = "Abundance")
  
  # Calculate mean abundance for each genus per stage
  grouped_df <- filtered_df_long %>%
    group_by(Genus, pathologic_stage_label) %>%
    summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")
  
  # Convert pathologic_stage_label to factor with ordered levels
  grouped_df$pathologic_stage_label <- factor(grouped_df$pathologic_stage_label, levels = c(1, 2, 3, 4))
  
  # Plot
  plot <- ggplot(grouped_df, aes(x = Genus, y = mean_abundance, fill = pathologic_stage_label)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Genus", y = "Mean Abundance", title = "Abundance of Genus' with Smallest P-Values from Kruskal Test") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set1", name = "Tumor Stage", labels = c("Stage 1", "Stage 2", "Stage 3", "Stage 4")) +  
    theme(legend.position = "top") + 
    ylim(0, NA)  # Ensure bars start at 0 level
  
  # Save plot to PDF
  pdf(output_file, width = 10, height = 6)  # Open PDF device
  print(plot)
  dev.off()  # Close PDF device
}

# Call the function with your input dataframe, min5_kruskal_taxa, and desired output file name
filter_and_plot_abundance(input_df = kraken_genusabundance, min5_kruskal_taxa = min5_kruskal_taxa, output_file = "min5_abundance_kruskal.pdf")
filter_and_plot_abundance(input_df = kraken_genusabundance, min5_kruskal_taxa = spear_maxabs_top5_taxa, output_file = "top5_spear_abundance_genus.pdf")

