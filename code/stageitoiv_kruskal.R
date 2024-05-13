kraken_metaCOADg_stageIIV <- kraken_meta_COAD_genus %>%
  filter(pathologic_stage_label %in% c("Stage I", "Stage IV"))

colnames(kraken_metaCOADg_stageIIV)[colnames(kraken_metaCOADg_stageIIV) == "...1"] <- "id"


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
  input_df <- input_df[, !grepl("pathologic_stage_label\\..*", colnames(input_df))]
  input_df <- input_df[, c("id", "pathologic_stage_label", setdiff(names(input_df), c("id", "pathologic_stage_label")))]
  input_df <- input_df[, -1]
  input_df <- input_df[, !grepl("pathologic_stage_label\\..*", colnames(input_df))]
  
  kruskal_results <- lapply(input_df[, -1], function(x) {
    kruskal_result <- kruskal.test(x ~ pathologic_stage_label, data = input_df)
    return(kruskal_result$p.value)
  })
  
  kruskal_results_df <- data.frame(taxon = names(kruskal_results), p_value = unlist(kruskal_results))
  
   #Apply FDR correction
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
 
  
  stage_colors <- c("Stage I" = "blue", "Stage IV" = "#FC4703")
  
  # Reorder the factor levels for pathologic_stage_label
  filtered_df_long$pathologic_stage_label <- factor(filtered_df_long$pathologic_stage_label, levels = c("Stage I", "Stage IV"))
  
  # Create the plot with manual color scale
  kruskal_plot <- ggplot(top_taxa_df, aes(x = reorder(Genus, match(Genus, kruskal_results_df$taxon)), y = Abundance, color = pathologic_stage_label)) +
    geom_boxplot(width = 0.75) +
    labs(x = "Genus", y = "Abundance", title = "Comparison of Taxon Abundance Across Stages for Taxon with Lowest P-value from Kruskal-Wallis Test") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 15, margin = margin(t = 30), color = "black"),
          axis.title.y = element_text(size = 15, margin = margin(r = 5), color="black"),
          axis.text.x = element_text(margin = margin(t=5), size = 15, color = "black"),  # Increase x-axis label size
          axis.text.y = element_text(margin = margin(r=5), size = 15, color= "black"),  # Increase y-axis label size
          axis.line = element_line(color = "black"),  # Change color of axis lines
          panel.grid = element_blank(),  # Remove all grid lines
          legend.text = element_text(size = 12),  # Increase legend text size
          legend.title = element_blank(),  # Increase legend title size
          legend.position = "right",  # Position legend at the top
          plot.title = element_text(size = 14.25, hjust = 0.29),
          legend.key.size = unit(2, "lines"),  # Increase space between legend elements
          plot.margin = margin(40, 20, 30, 30)) +  # Increase space below x-axis label and above graph
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +  
    scale_color_manual(values = stage_colors, name = "Tumor Stage", labels = names(stage_colors)) +
    guides(shape = guide_legend(override.aes = list(size = 40)))
  
  #Wilcos-rank-Sum test
  wilcox_results <- lapply(input_df[, -1], function(x) {
    wilcox_result <- wilcox.test(x ~ pathologic_stage_label, data = input_df)
    return(wilcox_result$p.value)
  })
  
  wilcox_results_df <- data.frame(taxon = names(wilcox_results), p_value_wilcox = unlist(wilcox_results))
  
  # Apply FDR correction for Wilcoxon Rank-Sum Test
  wilcox_results_df$p_adjust_wilcox <- p.adjust(wilcox_results_df$p_value_wilcox, method = "BH")
  wilcox_results_df <- wilcox_results_df %>%
    arrange(p_adjust_wilcox)
  print(wilcox_results_df)
  wilcox_results_df$rank_wilcox <- seq_len(nrow(wilcox_results_df))

  top_taxa_wilcox <- wilcox_results_df %>%
    slice_head(n = n_top) %>%
    pull(taxon)
  
  filtered_df_wilcox <- input_df %>%
    select(pathologic_stage_label, all_of(top_taxa_wilcox))
  
  filtered_df_long_w <- filtered_df_wilcox %>%
    pivot_longer(cols = -pathologic_stage_label, names_to = "Genus", values_to = "Abundance", values_drop_na = TRUE)
  
  # Filter the top taxa
  top_taxa_df_wilcox <- filtered_df_long_w %>%
    filter(Genus %in% top_taxa_wilcox)
  
  filtered_df_long_w$pathologic_stage_label <- factor(filtered_df_long_w$pathologic_stage_label, levels = c("Stage I", "Stage IV"))
  
  wilcox_plot <- ggplot(top_taxa_df_wilcox, aes(x = reorder(Genus, match(Genus, wilcox_results_df$taxon)), y = Abundance, color = pathologic_stage_label)) +
    geom_boxplot(width = 0.75) +
    labs(x = "Genus", y = "Abundance", title = "Comparison of Taxon Abundance Across Pathologic Stages for Taxon with Lowest P-value from Wilcoxon Rank-Sum Test") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 15, margin = margin(t = 5), color = "black"),
          axis.title.y = element_text(size = 15, margin = margin(r = 5), color="black"),
          axis.text.x = element_text(margin = margin(t=5), size = 15, color = "black"),  # Increase x-axis label size
          axis.text.y = element_text(margin = margin(r=5), size = 15, color= "black"),  # Increase y-axis label size
          axis.line = element_line(color = "black"),  # Change color of axis lines
          panel.grid = element_blank(),  # Remove all grid lines
          legend.text = element_text(size = 12),  # Increase legend text size
          legend.title = element_blank(),  # Increase legend title size
          legend.position = "right",  # Position legend at the top
          plot.title = element_text(size = 15, hjust = 0.2),
          legend.key.size = unit(2, "lines"),  # Increase space between legend elements
          plot.margin = margin(28, 20, 18, 18)) +  # Increase space below x-axis label and above graph
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +  
    scale_color_manual(values = stage_colors, name = "Tumor Stage", labels = names(stage_colors)) +
    guides(shape = guide_legend(override.aes = list(size = 40)))
  
  
  input_df$pathologic_stage_label <- as.numeric(factor(input_df$pathologic_stage_label, levels = c("Stage I","Stage IV")))
  
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
  
  top_taxa <- spearman_results_df %>%
    slice_head(n=2)
  spearman_results_df$correlation_sign <- ifelse(spearman_results_df$correlation > 0, "positive", "negative")
  print(spearman_results_df)
  # Volcano plot with labeled top ranked taxa
  volcano_plot_labeled <- ggplot(spearman_results_df, aes(x = correlation, y = -log10(p_adjust), color = correlation_sign)) +
    geom_point(alpha = 0.9) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_text(data = top_taxa, aes(label = taxon), hjust = 0.5, vjust = -0.7, size = 5, color = "black", bg = "transparent") +
    labs(title = expression(paste("Spearman's ", rho, " Correlation between Stage I and IV and Genus Abundance")),
         x = expression(paste("Spearman's ", rho)),
         y = "-log10(Adjusted p-value)",
         color = "Correlation") +
    scale_color_manual(values = c("positive" = "#EAB606", "negative" = "#006400"), 
                       labels = c("positive" = "Positive", "negative" = "Negative"))  +
    theme(
      plot.title = element_text(size = 15, margin = margin(b=5), hjust = 0.9),
      legend.position = "right",
      legend.title = element_text(size = 12), 
      legend.key.size = unit(1, "lines"),   # Adjust the size of legend keys
      legend.text = element_text(size = 12), # Adjust the size of legend text
      axis.title.y = element_text(size = 15, margin = margin(r = 2), color = "black"),  # Adjust the size of y-axis title
      axis.title.x = element_text(size = 15, margin = margin(t = 5), color = "black"),  # Adjust the size of x-axis title
      axis.line = element_line(color = "black"),
      plot.margin = margin(15, 10, 10, 10)) +  # Increase space below x-axis label and above graph
    xlim(c(-0.3, 0.3)) +  # Adjust the limits according to your data range
    theme(panel.background = element_blank())  +
    guides(color = guide_legend(override.aes = list(size = 4)))  # Adjust the size of legend keys (dots)
  # Remove plot background
  pdf(output_file, width = 5, height = 6)
  # Print the plot
  print(volcano_plot_labeled)
  print(kruskal_plot)
  print(wilcox_plot)
  
  dev.off()
  
  
  
  return(spearman_results_df)
}


wilcox_results <- analysis_stageI_stageIV(input_df = kraken_COAD_genus_stageIIV,
                        metadata_df = kraken_metaCOADg_stageIIV,            
                        n_top = 5,
                        output_file = "voomsnm_stageiiv_wilcoxkruskalspear.pdf")



library(dplyr)

norm_stageIIv_list <- lapply(data_list_with_stages, function(df) {
  df %>% filter(pathologic_stage_label %in% c("Stage I", "Stage IV"))
})

spear_comb_IIV <- list()

for (method in names(norm_stageIIv_list)){
  current_df <- norm_stageIIv_list[[method]]
  filtered_meta <- kraken_metaCOADg_stageIIV %>% 
    filter(id %in% current_df$id)
  current_df <- current_df %>% 
    select(-pathologic_stage_label)
  variable_name = paste0("spearstageIIv", method)
  output_file = paste0("spearstageIIv_", method, ".pdf")
  variable_name <- analysis_stageI_stageIV(input_df = current_df,
                          metadata_df = filtered_meta,
                          n_top = 5,
                          output_file = output_file)
  spear_comb_IIV[[method]] <- variable_name
  
}

combined_spearman_IIV <- do.call(rbind, spear_comb_IIV)
combined_spearman_IIV$Method <- sub("\\..*", "", rownames(combined_spearman_IIV))
combined_spearman_IIV <- combined_spearman_IIV[!combined_spearman_IIV$Method %in% c("RLE", "logcpm"), ]
print(unique(combined_spearman_IIV$Method))

ranksum_spearman_IIV <- combined_spearman_IIV %>%
  group_by(taxon) %>%
  filter(n_distinct(Method) == 10) %>%  # Filters taxon shared among all methods
  summarise(total_rank = sum(rank)) %>%
  arrange(total_rank)

top_taxa <- ranksum_spearman_IIV$taxon[1:10]  # Selecting top 10 taxa

# Filter the combined PCA dataset to include only the top taxa
subset_data_IIV <- combined_spearman_IIV %>%
  filter(taxon %in% top_taxa)

# Reorder the 'Method' factor variable
subset_data_IIV$Method <- factor(subset_data_IIV$Method, 
                             levels = c("Raw", "VoomSNM", "GMPR", "UQ", "TSS", 
                                        "RLEpos", "MED", "CLR", "CLRpos", "CSS"))

heatmap_plot <- ggplot(subset_data_IIV, aes(x = Method, y = factor(taxon, levels = rev(top_taxa)), fill = abs(correlation))) +
  geom_tile(color = "#ECECEC", linewidth = 0.2) +
  scale_fill_gradient(low = "white", high = "#36617B") +
  theme_minimal() +
  labs(title = expression(paste("Taxa with top Spearman's ", rho, " between Taxon Abundance and Stage I and IV Pooled Across Different Normalized Datasets")),
       x = "Normalization Method", y = "Genus", fill = "Correlation") +
  theme(axis.text.x = element_text(size = 7),  # Adjust the size of x-axis text
        axis.title.x = element_text(size = 8, margin = margin(t = 5), color = "black"),  # Adjust the size of x-axis title
        axis.title.y = element_text(size = 8, margin = margin(r = 2), color = "black"),  # Adjust the size of y-axis title
        axis.text.y = element_text(size = 7),  # Adjust the size of y-axis text
        legend.text = element_text(size = 8),  # Adjust the size of legend text
        legend.title = element_text(size = 8),  # Adjust the size of legend title
        legend.position = "right",  # Position legend
        plot.title = element_text(size = 9, margin = margin(b = 5), hjust = 0.5),  # Adjust the size of plot title
        legend.key.size = unit(1, "lines"))  # Adjust the space between legend elements

pdf("heatmap_spear_normtech_stageIIV.pdf", width = 7.5, height = 2)
print(heatmap_plot)
dev.off()


output_file <- "combined_volcano_plots.pdf"  # Define output file name

combined_plots <- list()  # Initialize a list to store combined plots

for (method in names(norm_stageIIv_list)) {
  if (method %in% c("GMPR", "UQ", "CLR", "TSS")) {
    current_df <- norm_stageIIv_list[[method]]
    filtered_meta <- kraken_metaCOADg_stageIIV %>% 
      filter(id %in% current_df$id)
    current_df <- current_df %>% 
      select(-pathologic_stage_label)
    
    # Preparing Data
    colnames(current_df)[1] <- "id"
    colnames(filtered_meta)[1] <- "id"
    input_df <- merge(current_df, filtered_meta[c('id','pathologic_stage_label')], by ='id',all.x=TRUE)
    if (any(grepl("pathologic_stage_label\\..", colnames(input_df)))) {
      colnames(input_df) <- make.unique(gsub("pathologic_stage_label\\..*", "pathologic_stage_label", colnames(input_df)))
    }
    input_df <- input_df[, !grepl("pathologic_stage_label\\..*", colnames(input_df))]
    input_df <- input_df[, c("id", "pathologic_stage_label", setdiff(names(input_df), c("id", "pathologic_stage_label")))]
    input_df <- input_df[, -1]
    input_df <- input_df[, !grepl("pathologic_stage_label\\..*", colnames(input_df))]
    
    input_df$pathologic_stage_label <- as.numeric(factor(input_df$pathologic_stage_label, levels = c("Stage I","Stage IV")))
    
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
      slice_head(n=2)
    
    spearman_results_df$correlation_sign <- ifelse(spearman_results_df$correlation > 0, "positive", "negative")
    
    # Volcano plot with labeled top ranked taxa
    volcano_plot_labeled <- ggplot(spearman_results_df, aes(x = correlation, y = -log10(p_adjust), color = correlation_sign)) +
      geom_point(alpha = 0.9) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      geom_text(data = top_taxa, aes(label = taxon), hjust = 0.5, vjust = -0.7, size = 5, color = "black", bg = "transparent") +
      labs(title = expression(paste("Spearman's ", rho, " Correlation between Stage I and IV and Genus Abundance")),
           x = expression(paste("Spearman's ", rho)),
           y = "-log10(Adjusted p-value)",
           color = "Correlation") +
      scale_color_manual(values = c("positive" = "#EAB606", "negative" = "#006400"), 
                         labels = c("positive" = "Positive", "negative" = "Negative"))  +
      theme(
        plot.title = element_text(size = 15, margin = margin(b=5), hjust = 0.9),
        legend.position = "right",
        legend.title = element_text(size = 12), 
        legend.key.size = unit(1, "lines"),   # Adjust the size of legend keys
        legend.text = element_text(size = 12), # Adjust the size of legend text
        axis.title.y = element_text(size = 15, margin = margin(r = 2), color = "black"),  # Adjust the size of y-axis title
        axis.title.x = element_text(size = 15, margin = margin(t = 5), color = "black"),  # Adjust the size of x-axis title
        axis.line = element_line(color = "black"),
        plot.margin = margin(15, 10, 10, 10)) +  # Increase space below x-axis label and above graph
      xlim(c(-0.3, 0.3)) +  # Adjust the limits according to your data range
      theme(panel.background = element_blank())  +
      guides(color = guide_legend(override.aes = list(size = 4)))  # Adjust the size of legend keys (dots)
    
    # Store plot in list
    combined_plots[[method]] <- volcano_plot_labeled
  }
}

# Combine plots into a single plot with facets
combined_plot <- do.call(grid.arrange, c(combined_plots, nrow = 2))

# Save the combined plot to a PDF file
output_file <- "combined_volcano_plots.pdf"
pdf(output_file, width = 12, height = 8)
print(combined_plot)
dev.off()

# Save the combined plot to a PDF file
output_file <- "combined_volcano_plots.pdf"
ggsave(output_file, combined_plot, width = 12, height = 8)







# Initialize a list to store combined plots
combined_plots <- list()

# Initialize variables to store minimum and maximum correlation and p-value across all plots
min_cor <- Inf
max_cor <- -Inf
min_p <- Inf
max_p <- -Inf

for (method in names(norm_stageIIv_list)) {
  if (method %in% c("GMPR", "UQ", "CLR", "TSS")) {
    current_df <- norm_stageIIv_list[[method]]
    filtered_meta <- kraken_metaCOADg_stageIIV %>% 
      filter(id %in% current_df$id)
    current_df <- current_df %>% 
      select(-pathologic_stage_label)
    
    # Preparing Data
    colnames(current_df)[1] <- "id"
    colnames(filtered_meta)[1] <- "id"
    input_df <- merge(current_df, filtered_meta[c('id','pathologic_stage_label')], by ='id',all.x=TRUE)
    if (any(grepl("pathologic_stage_label\\..", colnames(input_df)))) {
      colnames(input_df) <- make.unique(gsub("pathologic_stage_label\\..*", "pathologic_stage_label", colnames(input_df)))
    }
    input_df <- input_df[, !grepl("pathologic_stage_label\\..*", colnames(input_df))]
    input_df <- input_df[, c("id", "pathologic_stage_label", setdiff(names(input_df), c("id", "pathologic_stage_label")))]
    input_df <- input_df[, -1]
    input_df <- input_df[, !grepl("pathologic_stage_label\\..*", colnames(input_df))]
    
    input_df$pathologic_stage_label <- as.numeric(factor(input_df$pathologic_stage_label, levels = c("Stage I","Stage IV")))
    
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
      slice_head(n=2)
    
    spearman_results_df$correlation_sign <- ifelse(spearman_results_df$correlation > 0, "positive", "negative")
    
    # Volcano plot with labeled top ranked taxa
    volcano_plot_labeled <- ggplot(spearman_results_df, aes(x = correlation, y = -log10(p_adjust), color = correlation_sign)) +
      geom_point(alpha = 0.9) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      geom_text(data = top_taxa, aes(label = taxon), hjust = 0.5, vjust = -0.7, size = 5, color = "black", bg = "transparent") +
      labs(title = method,
           x = expression(paste("Spearman's ", rho)),
           y = "-log10(Adjusted p-value)",
           color = "Correlation") +
      scale_color_manual(values = c("positive" = "#EAB606", "negative" = "#006400"), 
                         labels = c("positive" = "Positive", "negative" = "Negative"))  +
      theme(
        plot.title = element_text(size = 15, margin = margin(b=5), hjust = 0.5),
        legend.position = "none",
        legend.title = element_text(size = 12), 
        legend.key.size = unit(1, "lines"),   # Adjust the size of legend keys
        legend.text = element_text(size = 12), # Adjust the size of legend text
        axis.title.y = element_text(size = 15, margin = margin(r = 2), color = "black"),  # Adjust the size of y-axis title
        axis.title.x = element_text(size = 15, margin = margin(t = 5), color = "black"),  # Adjust the size of x-axis title
        axis.line = element_line(color = "black"),
        plot.margin = margin(15, 10, 10, 10)) +  # Increase space below x-axis label and above graph
      xlim(c(-0.3, 0.3)) +  # Adjust the limits according to your data range
      theme(panel.background = element_blank())  +
      guides(color = guide_legend(override.aes = list(size = 4)))  # Adjust the size of legend keys (dots)
    
    # Store plot in list
    combined_plots[[method]] <- volcano_plot_labeled
    
    
    min_cor <- min(min_cor, min(spearman_results_df$correlation))
    max_cor <- max(max_cor, max(spearman_results_df$correlation))
    min_p <- min(min_p, min(-log10(spearman_results_df$p_adjust)))
    max_p <- max(max_p, max(-log10(spearman_results_df$p_adjust)))
    
    
  }
}

common_xlim <- c(min_cor, max_cor)
common_ylim <- c(min_p, max_p)

# Modify each plot to use common axis limits
for (method in names(combined_plots)) {
  combined_plots[[method]] <- combined_plots[[method]] +
    xlim(common_xlim) +
    ylim(common_ylim)
}

# Combine plots into a single plot with facets
combined_plot <- do.call(grid.arrange, c(combined_plots, ncol = length(combined_plots)))

# Save the combined plot to a PDF file
output_file <- "combined_volcano_plots.pdf"
ggsave(output_file, combined_plot, width = 20, height = 8)

