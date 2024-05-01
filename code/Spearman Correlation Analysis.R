perform_spearman_and_plot_abundance <- function(input_df, metadata_df, n_top, output_file) {
  # Preparing Data
  colnames(input_df)[1] <- "id"
  colnames(metadata_df)[1] <- "id"
  merged_df <- merge(input_df, metadata_df[c('id', 'pathologic_stage_label')], by = 'id', all.x = TRUE)
  filtered_df <- merged_df[, c("id", "pathologic_stage_label", setdiff(names(merged_df), c("id", "pathologic_stage_label")))]
  filtered_df <- filtered_df[, -1]
  filtered_df$pathologic_stage_label <- as.numeric(factor(filtered_df$pathologic_stage_label, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")))
  
  # Calculate Spearman correlation for each taxon
  spearman_correlation <- apply(filtered_df[, -c(1, ncol(filtered_df))], 2, function(x) {
    cor(filtered_df$pathologic_stage_label, x, method = "spearman")
  })
  
  # Convert results to dataframe
  spearman_results_df <- data.frame(taxon = names(spearman_correlation), correlation = unlist(spearman_correlation))
  
  # Apply FDR correction
  spearman_results_df$p_adjust <- p.adjust(abs(spearman_results_df$correlation), method = "fdr")
  
  # Arrange by absolute correlation value
  spearman_results_df <- spearman_results_df %>%
    arrange(-abs(correlation))
  
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
  
  # Plot
  plot <- ggplot(filtered_df_long, aes(x = reorder(Genus, match(Genus, spearman_results_df$taxon)), y = Abundance, color = pathologic_stage_label)) +
    geom_jitter(position = position_dodge(width = 0.75)) +
    labs(x = "Genus", y = "Abundance", title = paste("Abundance of Top", n_top, "Genus with Highest Absolute Spearman Correlation with Tumor Stage")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_brewer(palette = "Set1", name = "Tumor Stage") +
    theme(legend.position = "top")
  
  # Save plot to PDF
  pdf(output_file, width = 12, height = 8)
  print(plot)
  dev.off()  # Close PDF device
  
  # Print the list of top taxa and their corresponding Spearman coefficients
  cat("List of Top 25 Taxa with Highest Absolute Spearman Correlation:\n")
  top_taxa_correlation <- spearman_results_df %>%
    slice_head(n = 25)
  
  # Return the top taxa with highest absolute correlation
  return(top_taxa)
}




perform_spearman_and_plot_abundance(input_df = kraken_COAD_genus,
                                    metadata_df = kraken_meta_COAD_genus,
                                    n_top = 5,
                                    output_file = "totalkraken_voomsnm_absmax5_spearman.pdf")

perform_spearman_and_plot_abundance(input_df = kraken_COAD_IlGA_UNC_g,
                                    metadata_df = kraken_metaCOAD_IlGA_UNC_g,
                                    n_top = 5,
                                    output_file = "uncrna_voomsnm_absmax5_spearman.pdf")




# Call the function with your input dataframe, number of top taxa, and desired output file name
perform_spearman_and_plot_abundance(input_df = kraken_COAD_genus,
                                    metadata_df = kraken_meta_COAD_genus,            
                                    n_top = 5,
                                    output_file = "totalkraken_voomsnm_absmax5_spearman.pdf")

perform_spearman_and_plot_abundance(input_df = `kraken_COADg_UniversityofNorthCarolina_RNA-Seq`,
                                                              metadata_df = `kraken_meta_UniversityofNorthCarolina_RNA-Seq`,                 
                                                              n_top = 5,
                                                              output_file = "unc_voomsnm_absmax5_spearman.pdf")

perform_spearman_and_plot_abundance(input_df = `kraken_COADg_BaylorCollegeofMedicine_WGS`,
                                                                  metadata_df = `kraken_meta_BaylorCollegeofMedicine_WGS`,                 
                                                                  n_top = 5,
                                                                  output_file = "baylor_voomsnm_absmax5_spearman.pdf")

perform_spearman_and_plot_abundance(input_df = `kraken_COADg_HarvardMedicalSchool_WGS`,
                                                                   metadata_df = `kraken_meta_HarvardMedicalSchool_WGS`,                 
                                                                   n_top = 5,
                                                                   output_file = "harvard_voomsnm_absmax5_spearman.pdf")


#Looking at Some based on Tissue Source
perform_spearman_and_plot_abundance(input_df = `kraken_COADg_ChristianaHealthcare`,
                                                                            metadata_df = `kraken_meta_ChristianaHealthcare`,                 
                                                                            n_top = 5,
                                                                            output_file = "christiana_voomsnm_absmax5_spearman.pdf")

perform_spearman_and_plot_abundance(input_df = `kraken_COADg_MSKCC`,
                                                                       metadata_df = `kraken_meta_MSKCC`,                 
                                                                       n_top = 5,
                                                                       output_file = "MSKCC_voomsnm_absmax5_spearman.pdf")


