#Kruskal + Spear Analysis for Submitting Center
colnames(kraken_COAD_genus)[1] <- "id"
colnames(kraken_metaCOAD)[1] <- "id"
split_metadata_submittingcenter <- split(kraken_metaCOAD, kraken_metaCOAD$data_submitting_center_label)
unique_centers <- unique(kraken_metaCOAD$data_submitting_center_label)

library(tidyr)


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


library(dplyr)

# Combine datasets and add a dataset column
comsubmit_spear <- bind_rows(
  `spear_Harvard Medical School` %>% mutate(dataset = "Harvard Medical School"),
  `spear_University of North Carolina` %>% mutate(dataset = "University of North Carolina")
)

# Calculate total rank for each taxa
ranksum_spear_submit <- comsubmit_spear %>%
  group_by(taxon) %>%
  summarise(total_rank = sum(rank))
ranksum_spear_submit <- ranksum_spear_submit %>%
  arrange(total_rank)  # Arrange in descending order of total rank

# Select top N taxa (adjust N as needed)
top_taxa <- ranksum_spear_submit$taxon[1:5]  # Selecting top 10 taxa

# Filter the combined PCA dataset to include only the top taxa
subset_data <- comsubmit_spear %>%
  filter(taxon %in% top_taxa)



# Create the heatmap plot with adjusted saturation
heatmap_plot <- ggplot(subset_data, aes(x = factor(taxon, levels = top_taxa), y = dataset, fill = p_adjust)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "white", high = "#0D3601", 
                      limits = c(0.1, 0.3),  # Set the range of values
                      breaks = seq(0.1, 0.5, 0.05),  # Define breaks for color scale
                      trans = "sqrt") +  # Adjust saturation with transformation function
  theme_minimal() +
  labs(title = "Top Ranking Taxa Across Combined PCA Datasets",
       x = "Taxa", y = "Dataset", fill = "Correlation") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14))

# Save the plot as a PDF with larger dimensions
pdf("heatmap_plot.pdf", width = 20, height = 10)
print(heatmap_plot)
dev.off()


#BASED ON PCA ONLY LOOKING AT UNC AND HARVARD SUBMITTING
comsubmit_kruskal <- rbind(`kruskal_Harvard Medical School`, `kruskal_University of North Carolina`)
ranksum_kruskal_submit <- comsubmit_kruskal %>%
  group_by(taxon) %>%
  summarise(total_rank = sum(rank))
ranksum_kruskal_submit <- ranksum_kruskal_submit %>%
  arrange(total_rank)



