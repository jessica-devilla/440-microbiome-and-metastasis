library(tidyr)

perform_spearman_and_plot_abundance <- function(input_df, metadata_df, n_top, output_file) {
  # Preparing Data
  colnames(input_df)[1] <- "id"
  colnames(metadata_df)[1] <- "id"
  merged_df <- merge(input_df, metadata_df[c('id', 'pathologic_stage_label')], by = 'id', all.x = TRUE)
  if (any(grepl("pathologic_stage_label\\..", colnames(merged_df)))) {
    colnames(merged_df) <- gsub("pathologic_stage_label\\..*", "pathologic_stage_label", colnames(merged_df))
  }
  merged_df <- merged_df[, !grepl("pathologic_stage_label\\..*", colnames(merged_df))]
  filtered_df <- merged_df[, c("id", "pathologic_stage_label", setdiff(names(merged_df), c("id", "pathologic_stage_label")))]
  filtered_df <- filtered_df[, -1]
  filtered_df$pathologic_stage_label <- as.numeric(factor(filtered_df$pathologic_stage_label, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")))
  
  # Calculate Spearman correlation for each taxon
  spearman_correlation <- apply(filtered_df[, -c(1, ncol(filtered_df))], 2, function(x) {
    cor(filtered_df$pathologic_stage_label, x, method = "spearman")
  })
  
  # Convert results to dataframe
  spearman_results_df <- data.frame(taxon = names(spearman_correlation), correlation = unlist(spearman_correlation))
  
  spearman_results_df <- spearman_results_df %>%
    arrange(-abs(correlation))
  
  spearman_results_df$rank <- seq_len(nrow(spearman_results_df))
  
  # Calculate Spearman correlation and p-values for each taxon
  spearman_results <- lapply(filtered_df[, -c(1, ncol(filtered_df))], function(x) {
    cor_test_result <- cor.test(filtered_df$pathologic_stage_label, x, method = "spearman", exact = FALSE)
    return(list(correlation = cor_test_result$estimate, p_value = cor_test_result$p.value))
  })
  
  # Create a dataframe from the results
  spearman_results_df <- data.frame(
    taxon = names(spearman_results),
    correlation = sapply(spearman_results, function(x) x$correlation),
    p_value = sapply(spearman_results, function(x) x$p_value)
  )
  
  # Rank the taxa by their absolute correlation coefficient
  spearman_results_df <- spearman_results_df %>%
    mutate(rank = rank(-abs(correlation)))
  
  # Adjust p-values using the false discovery rate (FDR) method
  spearman_results_df$p_adjust <- p.adjust(spearman_results_df$p_value, method = "fdr")
  
  # Sort the dataframe by rank
  spearman_results_df <- spearman_results_df %>%
    arrange(rank)
  
  
  # Select top n taxa with highest absolute correlation after correction
  top_taxa <- spearman_results_df %>%
    slice_head(n = n_top) %>%
    pull(taxon)
  
  # Filter the input dataframe to include only the top taxa
  filtered_df <- merged_df %>%
    select(pathologic_stage_label, all_of(top_taxa))
  
  # Reshape the dataframe to long format for easy plotting
  filtered_df_long <- filtered_df %>%
    pivot_longer(cols = -pathologic_stage_label, names_to = "Genus", values_to = "Abundance")
  
  stage_colors <- c("Stage I" = "blue", "Stage II" = "#8070FE", "Stage III" = "#EAB606", "Stage IV" = "#FC4703")
  
  pdf(output_file, width = 20, height = 10)
  # Plot
  plot <- ggplot(filtered_df_long, aes(x = reorder(Genus, match(Genus, spearman_results_df$taxon)), y = Abundance, color = pathologic_stage_label)) +
    geom_boxplot(width = 0.75) + 
    labs(x = "Genus", y = "Abundance", title = "Genus Abundance from Top 5 Taxa in Spearman Correlation") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 22, margin = margin(t = 20), color = "black"),
          axis.title.y = element_text(size = 22, margin = margin(r = 15), color="black"),
          axis.text.x = element_text(margin = margin(t=15), size = 20, color = "black"),  # Increase x-axis label size
          axis.text.y = element_text(margin = margin(r=15), size = 20, color= "black"),  # Increase y-axis label size
          axis.line = element_line(color = "black"),  # Change color of axis lines
          panel.grid = element_blank(),  # Remove all grid lines
          legend.text = element_text(size = 18),  # Increase legend text size
          legend.title = element_blank(),  # Increase legend title size
          legend.position = "right",  # Position
          plot.title = element_text(size = 30, margin = margin(b=20), hjust = 0.5),
          legend.key.size = unit(2, "lines"),  # Increase space between legend elements
          plot.margin = margin(30, 20, 20, 20)) +  # Increase space below x-axis label and above graph
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +  
    scale_color_manual(values = stage_colors, name = "Tumor Stage", labels = names(stage_colors)) +
    guides(shape = guide_legend(override.aes = list(size = 40)))
  
  # Save plot to PDF
  print(plot)
  dev.off()  # Close PDF device
  
  filtered_df$pathologic_stage_label <- as.numeric(factor(filtered_df$pathologic_stage_label, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")))
  
  pearson_correlation <- apply(filtered_df[, -c(1, ncol(filtered_df))], 2, function(x) {
    cor(filtered_df$pathologic_stage_label, x, method = "pearson")
  })
  
  # Convert results to dataframe
  pearson_results_df <- data.frame(taxon = names(pearson_correlation), pearson_correlation = unlist(pearson_correlation))
  
  # Calculate p-values for Pearson correlation
  pearson_p_values <- apply(filtered_df[, -c(1, ncol(filtered_df))], 2, function(x) {
    cor_test_result <- cor.test(filtered_df$pathologic_stage_label, x, method = "pearson", exact = FALSE)
    return(cor_test_result$p.value)
  })
  
  # Convert p-values to dataframe
  pearson_p_values_df <- data.frame(taxon = names(pearson_p_values), pearson_p_value = unlist(pearson_p_values))
  
  # Adjust p-values using the false discovery rate (FDR) method
  pearson_p_values_df$pearson_p_adjust <- p.adjust(pearson_p_values_df$pearson_p_value, method = "fdr")
  print(pearson_p_values_df)
  # Rank the taxa by their absolute correlation coefficient
  pearson_results_df <- pearson_results_df %>%
    mutate(rank = rank(-abs(pearson_correlation)))
  
  print(pearson_results_df)
}

# Call the function with your input dataframe, number of top taxa, and desired output file name
spear_vnm_total <- perform_spearman_and_plot_abundance(input_df = kraken_COAD_genus,
                                                       metadata_df = kraken_meta_COAD_genus,            
                                                       n_top = 5,
                                                       output_file = "totalkraken_voomsnm_absmax5_spearman.pdf")

print(max(kraken_COAD_genus$Oscillatoria))



library(ggplot2)
pdf("volcano_spearman.pdf")
top_taxa <- spear_vnm_total %>%
  slice_head(n = 1)  # Change the number to select the top N taxa

# Convert sign of correlation to a factor with two levels
spear_vnm_total$correlation_sign <- ifelse(spear_vnm_total$correlation > 0, "positive", "negative")

# Volcano plot with labeled top ranked taxa
volcano_plot_labeled <- ggplot(spear_vnm_total, aes(x = correlation, y = -log10(p_adjust), color = correlation_sign)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text(data = top_taxa, aes(label = taxon), hjust = -0.2, vjust = -0.5, size = 4, color = "black", bg = "transparent") +
  labs(title = expression(paste("Spearman's ", rho, " Correlation between Stage and Genus within Colon Adenocarcinoma Tumor Microbiome")),
       x = expression(paste("Spearman's ", rho)),
       y = "-log10(Adjusted p-value)",
       color = "Correlation") +
  scale_color_manual(values = c("positive" = "#EAB606", "negative" = "#8070FE"), 
                     labels = c("positive" = "Positive", "negative" = "Negative")) +
  theme(
    plot.title = element_text(size = 10, margin = margin(b=10), hjust = 0.5),
    legend.position = "bottom",
    axis.line = element_line(color = "black")) +  # Add axis lines for both x and y axes
  xlim(c(-0.3, 0.3)) +  # Adjust the limits according to your data range
  theme(panel.background = element_blank())  # Remove plot background

# Print the plot
print(volcano_plot_labeled)

dev.off()

