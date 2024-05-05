kraken_metaCOADg_stageIIV <- kraken_meta_COAD_genus %>%
  filter(pathologic_stage_label %in% c("Stage I", "Stage IV"))

kraken_COAD_genus_stageIIV <- kraken_COAD_genus %>%
  filter(...1%in%kraken_metaCOADg_stageIIV$id)


analysis_stageI_stageIV <- function(input_df, metadata_df, n_top, output_file) {
  #Preparing Data
  colnames(input_df)[1] <- "id"
  colnames(metadata_df)[1] <- "id"
  input_df <- merge(input_df, metadata_df[c('id','pathologic_stage_label')], by ='id',all.x=TRUE)
  if (any(grepl("pathologic_stage_label\\..", colnames(input_df)))) {
    colnames(input_df) <- make.unique(gsub("pathologic_stage_label\\..*", "pathologic_stage_label", colnames(input_df)))
  }
  input_df <- input_df[, c("id", "pathologic_stage_label", setdiff(names(input_df), c("id", "pathologic_stage_label")))]
  input_df <- input_df[, -1]
  input_df <- input_df[, !grepl("pathologic_stage_label\\..*", colnames(input_df))]
  
  kruskal_results <- lapply(input_df[, -1], function(x) {
    kruskal_result <- kruskal.test(x ~ pathologic_stage_label, data = input_df)
    return(kruskal_result$p.value)
  })
  
  kruskal_results_df <- data.frame(taxon = names(kruskal_results), p_value = unlist(kruskal_results))
  
  # Apply FDR correction
  kruskal_results_df$p_adjust <- p.adjust(kruskal_results_df$p_value, method = "BH")
  kruskal_results_df <- kruskal_results_df %>%
    arrange(p_adjust)
  kruskal_results_df$rank <- seq_len(nrow(kruskal_results_df))
  
  kruskal_results_df$rank <- seq_len(nrow(kruskal_results_df))
  
  print(kruskal_results_df)
  # Create a separate dataframe for ranked taxa
  ranked_taxa_df <- kruskal_results_df %>%
    select(taxon, rank)
  
  # Select top n taxa with smallest p-values after correction
  top_taxa <- kruskal_results_df %>%
    slice_head(n = n_top) %>%
    pull(taxon)
  
  # Filter the input dataframe to include only the top taxa
  filtered_df <- input_df %>%
    select(pathologic_stage_label, all_of(top_taxa))
  
  filtered_df_long <- filtered_df %>%
    pivot_longer(cols = -pathologic_stage_label, names_to = "Genus", values_to = "Abundance", values_drop_na = TRUE)
  
  # Filter the top taxa
  top_taxa_df <- filtered_df_long %>%
    filter(Genus %in% top_taxa)
  print(top_taxa)
  
  stage_colors <- c("Stage I" = "blue", "Stage IV" = "#FC4703")
  
  # Reorder the factor levels for pathologic_stage_label
  filtered_df_long$pathologic_stage_label <- factor(filtered_df_long$pathologic_stage_label, levels = c("Stage I", "Stage IV"))
  
  # Create the plot with manual color scale
  kruskal_plot <- ggplot(top_taxa_df, aes(x = reorder(Genus, match(Genus, kruskal_results_df$taxon)), y = Abundance, color = pathologic_stage_label)) +
    geom_boxplot(width = 0.75) +
    labs(x = "Genus", y = "Abundance", title = "Genus Abunance Across Stages for Genus' with Greatest Variance between Stages") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 20, margin = margin(t = 10), color = "black"),
          axis.title.y = element_text(size = 20, margin = margin(r = 10), color="black"),
          axis.text.x = element_text(margin = margin(t=15), size = 15, color = "black"),  # Increase x-axis label size
          axis.text.y = element_text(margin = margin(r=15), size = 15, color= "black"),  # Increase y-axis label size
          axis.line = element_line(color = "black"),  # Change color of axis lines
          panel.grid = element_blank(),  # Remove all grid lines
          legend.text = element_text(size = 18),  # Increase legend text size
          legend.title = element_blank(),  # Increase legend title size
          legend.position = "right",  # Position legend at the top
          plot.title = element_text(size = 22, hjust = 0.5),
          legend.key.size = unit(2, "lines"),  # Increase space between legend elements
          plot.margin = margin(40, 20, 20, 20)) +  # Increase space below x-axis label and above graph
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +  
    scale_color_manual(values = stage_colors, name = "Tumor Stage", labels = names(stage_colors)) +
    guides(shape = guide_legend(override.aes = list(size = 40)))
  
  
  
  input_df$pathologic_stage_label <- as.numeric(factor(filtered_df$pathologic_stage_label, levels = c("Stage I","Stage IV")))
  
  # Calculate Spearman correlation for each taxon
  spearman_correlation <- apply(filtered_df[, -c(1, ncol(filtered_df))], 2, function(x) {
    cor(input_df$pathologic_stage_label, x, method = "spearman")
  })
  
  # Convert results to dataframe
  spearman_results_df <- data.frame(taxon = names(spearman_correlation), correlation = unlist(spearman_correlation))
  
  spearman_results_df <- spearman_results_df %>%
    arrange(-abs(correlation))
  
  spearman_results_df$rank <- seq_len(nrow(spearman_results_df))
  
  # Calculate Spearman correlation and p-values for each taxon
  spearman_results <- lapply(input_df[, -c(1, ncol(input_df))], function(x) {
    cor_test_result <- cor.test(input_df$pathologic_stage_label, x, method = "spearman", exact = FALSE)
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
  
  top_taxa <- spearman_results_df %>%
    slice_head(n = n_top)  # Change the number to select the top N taxa
  
  print(top_taxa)
  
  top_taxa <- spearman_results_df %>%
    slice_head(n=2)
  
  current_df$correlation_sign <- ifelse(current_df$correlation >= 0, "positive", "negative")
  
  # Volcano plot with labeled top ranked taxa
  volcano_plot_labeled <- ggplot(spearman_results_df, aes(x = correlation, y = -log10(p_adjust), color = correlation_sign)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_text(data = top_taxa, aes(label = taxon), hjust = -0.2, vjust = -0.5, size = 4, color = "black", bg = "transparent") +
    labs(title = bquote(paste("Spearman's ", rho, " between Stage and Genus in ", .(method), " on Tumor Microbiome Data")),
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
  
  pdf(output_file, width = 20, height = 10)
  print(kruskal_plot)
  print(volcano_plot_labeled)
  dev.off()  # Close PDF device
}


analysis_stageI_stageIV(input_df = kraken_COAD_genus_stageIIV,
                        metadata_df = kraken_metaCOADg_stageIIV,            
                        n_top = 5,
                        output_file = "voomsnm_stageiiv_kruskalspear.pdf")



library(dplyr)

norm_stageIIv_list <- lapply(data_list_with_stages, function(df) {
  df %>% filter(pathologic_stage_label %in% c("Stage I", "Stage IV"))
})

for (method in names(norm_stageIIv_list)){
  current_df <- norm_stageIIv_list[[method]]
  current_df$id <- rownames(current_df)
  filtered_meta <- kraken_metaCOADg_stageIIV %>% filter(id %in% current_df$id)
  current_df <- current_df %>% 
    select(-pathologic_stage_label)
  output_file = paste0("spearstageIIv_", method, ".pdf")
  analysis_stageI_stageIV(input_df = current_df,
                          metadata_df = filtered_meta,
                          n_top = 5,
                          output_file = output_file)
}

  