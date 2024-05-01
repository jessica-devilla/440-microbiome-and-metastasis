library(tidyr)


#Prep Voom-SNM Kraken Data
kraken_meta_COAD_genus <- kraken_metaCOAD

#Prep Voom-SNM Kraken UNC/IlluminaGA/RNA-seq Data
kraken_COAD_IlGA_UNC_g <- kraken_COADdata_RNASeq_IlluminaGA_UNC
colnames(kraken_COAD_IlGA_UNC_g) <- sub(".*__(.*)$", "\\1", colnames(kraken_COAD_IlGA_UNC_g))
kraken_metaCOAD_IlGA_UNC_g <- kraken_metaCOAD_RNASeq_IlluminaGA_UNC


#Function Preforms Kruskal Analysis, Plots the Minimum p-values of cohort according to n_top, and Plots Mean Abundance of Taxa w/ Min P-value
perform_kruskal_and_plot_abundance <- function(input_df, metadata_df, n_top, output_file) {
  #Preparing Data
  colnames(input_df)[1] <- "id"
  colnames(metadata_df)[1] <- "id"
  input_df <- merge(input_df, metadata_df[c('id','pathologic_stage_label')], by='id',all.x=TRUE)
  input_df <- input_df[, c("id", "pathologic_stage_label", setdiff(names(input_df), c("id", "pathologic_stage_label")))]
  input_df <- input_df[, -1]
  
  
  # Perform Kruskal-Wallis test on each taxon
  kruskal_results <- lapply(input_df[, -1], function(x) {
    kruskal_result <- kruskal.test(x ~ pathologic_stage_label, data = input_df)
    return(kruskal_result$p.value)
  })
  
  # Convert results to dataframe
  kruskal_results_df <- data.frame(taxon = names(kruskal_results), p_value = unlist(kruskal_results))

  # Apply FDR correction
  kruskal_results_df$p_adjust <- p.adjust(kruskal_results_df$p_value, method = "BH")
  kruskal_results_df <- kruskal_results_df %>%
    arrange(p_adjust)
  
  #kruskal_results_df$taxon <- gsub("_", " ", kruskal_results_df$taxon)
  
  kruskal_results_df$rank <- seq_len(nrow(kruskal_results_df))
  
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
    pivot_longer(cols = -pathologic_stage_label, names_to = "Genus", values_to = "Abundance")
  
  # Filter the top taxa
  top_taxa_df <- filtered_df_long %>%
    filter(Genus %in% top_taxa)
  
  stage_colors <- c("Stage I" = "red", "Stage II" = "blue", "Stage III" = "green", "Stage IV" = "purple")
  
  # Reorder the factor levels for pathologic_stage_label
  filtered_df_long$pathologic_stage_label <- factor(filtered_df_long$pathologic_stage_label, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
  
  # Create the plot with manual color scale
  plot <- ggplot(top_taxa_df, aes(x = reorder(Genus, match(Genus, kruskal_results_df$taxon)), y = Abundance, color = pathologic_stage_label)) +
    geom_jitter(position = position_dodge(width = 0.75)) +
    labs(x = "Genus", y = "Abundance", title = "Abundance of Genus' with Smallest P-Values from Kruskal Test") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +  # Replace "_" with " " in x-axis labels
    scale_color_manual(values = stage_colors, name = "Tumor Stage", labels = names(stage_colors)) +  
    theme(legend.position = "top")
  
  kruskal_results_df$taxon <- gsub("_", " ", kruskal_results_df$taxon)
  
  # Plot p-values
  plot_p_values <- ggplot(kruskal_results_df[1:n_top, ], aes(x = reorder(taxon, p_adjust), y = p_adjust)) +
    geom_point(stat = "identity", size = 3) +
    labs(x = "Taxon", y = "Adjusted p-value", title = paste("Minimum Kruskal-Wallis P-values (FDR corrected)")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(aes(label = round(p_adjust, 6)), vjust = -0.5) +
    ylim(0, NA)  # Ensure y-axis starts at 0
  
  # Save plots to PDF
  pdf(output_file, width = 12, height = 8)
  print(plot)
  print(plot_p_values)
  dev.off()  # Close PDF device
  
  # Return the top taxa with smallest p-values
  return(ranked_taxa_df)
}


# Call the function with your input dataframe, number of top taxa, and desired output file name
min5_kruskal_voomsnm <- perform_kruskal_and_plot_abundance(input_df = kraken_COAD_genus,
                                   metadata_df = kraken_meta_COAD_genus,            
                                   n_top = 5,
                                   output_file = "totalkraken_voomsnm_min5_kruskal.pdf")


