#Kruskal+ Spear Analysis for Tissue Source Site
colnames(kraken_COAD_genus)[1] <- "id"
colnames(kraken_metaCOAD)[1] <- "id"
split_meta_source <- split(kraken_metaCOAD, kraken_metaCOAD$tissue_source_site_label)
unique_source <- unique(kraken_metaCOAD$tissue_source_site_label)


for (source in unique_source) {
  dataset_name <- paste0("kraken_meta_", gsub(" ", "", source)) # Generate a unique variable name
  assign(dataset_name, split_meta_source[[source]]) # Assign the dataset to the variable
  
  current_dataset <- get(dataset_name)
  
  # Check if the number of samples is greater than 5
  if (nrow(current_dataset) > 10) {
    data_kraken_name <- paste0("kraken_COADg_", gsub(" ", "", source)) # Generate a unique variable name for kraken_COAD
    filtered_kraken <- kraken_COAD_genus %>%
      filter(id %in% current_dataset$id) # Replace ...1 with the appropriate column name from kraken_COAD
    
    assign(data_kraken_name, filtered_kraken) # Assign the filtered kraken_COAD dataset
    
    kruskal_name <- paste0("kruskaltsource_", gsub("","",source))
    output_file <- paste0(kruskal_name, ".pdf")
    
    kruskal <- perform_kruskal_and_plot_abundance(input_df = filtered_kraken,
                                                  metadata_df = split_meta_source[[source]],
                                                  n_top = 5,
                                                  output_file = output_file)
    assign(kruskal_name, kruskal)
    
    spear_name <- paste0("speartsource_", gsub("","",source))
    output_file2 <- paste0(spear_name, ".pdf")
    
    spear <- perform_spearman_and_plot_abundance(input_df = filtered_kraken,
                                                 metadata_df = split_meta_source[[source]],
                                                 n_top = 5,
                                                 output_file = output_file2)
    
    assign(spear_name, spear)
  } else {
    cat(paste("Skipping analysis for source", source, "due to insufficient samples\n"))
  }
}
