#Kruskal + Spear Analysis for Submitting Center
colnames(kraken_COAD_genus)[1] <- "id"
colnames(kraken_metaCOAD)[1] <- "id"
split_metadata_submittingcenter <- split(kraken_metaCOAD, kraken_metaCOAD$data_submitting_center_label)
unique_centers <- unique(kraken_metaCOAD$data_submitting_center_label)

for (center in unique_centers) {
  dataset_name <- paste0("kraken_meta_", gsub(" ", "", center)) # Generate a unique variable name
  split_metadata_submittingcenter <- split(kraken_metaCOAD, kraken_metaCOAD$data_submitting_center_label)
  assign(dataset_name, split_metadata_submittingcenter[[center]]) # Assign the dataset to the variable

  current_dataset <- get(dataset_name)
  
  data_kraken_name <- paste0("kraken_COADg_", gsub(" ", "", center)) # Generate a unique variable name for kraken_COAD
  filtered_kraken <- kraken_COAD_genus %>%
    filter(id %in% current_dataset$id) # Replace ...1 with the appropriate column name from kraken_COAD
  
  assign(data_kraken_name, filtered_kraken) # Assign the filtered kraken_COAD dataset
  kruskal_name <- paste0("kruskal_", gsub("","",center))
  output_file <- paste0(kruskal_name, ".pdf")

  
  kruskal <- perform_kruskal_and_plot_abundance(input_df = filtered_kraken,
                                                metadata_df = split_metadata_submittingcenter[[center]],
                                                n_top = 5,
                                                output_file = output_file)
  assign(kruskal_name, kruskal)
  
  
  spear_name <- paste0("spear_", gsub("","",center))
  output_file2 <- paste0(spear_name, ".pdf")
  
  spear <- perform_spearman_and_plot_abundance(input_df = filtered_kraken,
                                               metadata_df = split_metadata_submittingcenter[[center]],
                                               n_top = 5,
                                               output_file = output_file2)
  
  assign(spear_name, spear)
  
}


#BASED ON PCA ONLY COMBINING UNC AND HARVARD DATA
comsubmit_spear <- rbind(`spear_Harvard Medical School`, `spear_University of North Carolina`)
ranksum_spear_submit <- comsubmit_spear  %>%
  group_by(taxon) %>%
  summarise(total_rank = sum(rank))
ranksum_spear_submit <- ranksum_spear_submit %>%
  arrange(total_rank)


#BASED ON PCA ONLY LOOKING AT UNC AND HARVARD SUBMITTING
comsubmit_kruskal <- rbind(`kruskal_Harvard Medical School`, `kruskal_University of North Carolina`)
ranksum_kruskal_submit <- comsubmit_kruskal %>%
  group_by(taxon) %>%
  summarise(total_rank = sum(rank))
ranksum_kruskal_submit <- ranksum_kruskal_submit %>%
  arrange(total_rank)



