colon_metadata <- subset(kraken_metaCOAD, sample_type == "Primary Tumor")
Kraken_TCGA_Voom_SNM_Plate_Center_Filtering_Data <- kraken_COAD
install.packages("dplyr")  # Install the dplyr package
library(dplyr)

stageI_metadata <- subset(colon_metadata, pathologic_stage_label %in% c("Stage IA", "Stage IB","Stage I"))
stageII_metadata <- subset(colon_metadata, pathologic_stage_label %in% c("Stage IIA", "Stage IIB", "Stage II"))
stageIII_metadata <- subset(colon_metadata, pathologic_stage_label %in% c("Stage IIIA", "Stage IIIB", "Stage III"))
stageIV_metadata <- subset(colon_metadata, pathologic_stage_label %in% c("Stage IVA", "Stage IVB", "Stage IV"))

stageI_data <- Kraken_TCGA_Voom_SNM_Plate_Center_Filtering_Data %>%
  filter(`...1` %in% stageI_metadata$'...1')
stageII_data <- Kraken_TCGA_Voom_SNM_Plate_Center_Filtering_Data %>%
  filter(`...1` %in% stageII_metadata$'...1')
stageIII_data <- Kraken_TCGA_Voom_SNM_Plate_Center_Filtering_Data %>%
  filter(`...1` %in% stageIII_metadata$'...1')
stageIV_data <- Kraken_TCGA_Voom_SNM_Plate_Center_Filtering_Data %>%
  filter(`...1` %in% stageIV_metadata$'...1')


extract_unique_genus_names <- function(stage_data) {
  genus_names <- gsub("^.*g__([^\\.]+).*", "\\1", names(stage_data)[-1])
  genus_names <- genus_names[!grepl("^k__", genus_names)]
  genus_names <- genus_names[!grepl("^contaminant", genus_names)]
  unique_genus_names <- unique(genus_names)
  return(unique_genus_names)
}

# Extract unique genus names for each stage
unique_genus_names_stageI <- extract_unique_genus_names(stageI_data)
unique_genus_names_stageII <- extract_unique_genus_names(stageII_data)
unique_genus_names_stageIII <- extract_unique_genus_names(stageIII_data)
unique_genus_names_stageIV <- extract_unique_genus_names(stageIV_data)

calculate_genus_quantity <- function(stage_data, unique_genus_names, total_samples) {
  genus_quantity_list <- list()
  
  # Iterate over each unique genus name
  for (genus_name in unique_genus_names) {
    # Filter columns corresponding to the current genus
    genus_columns <- grep(paste0("\\.g__", genus_name), names(stage_data)[-1], value = TRUE, ignore.case = TRUE)
    
    # Calculate the total amount for the current genus across all samples
    total_amount_genus <- sum(stage_data[, genus_columns])
    
    # Normalize the total amount by the total number of samples
    normalized_quantity <- total_amount_genus / total_samples
    
    # Store the genus name and quantity in a list
    genus_quantity_list[[genus_name]] <- normalized_quantity
  }
  
  # Convert the list to a data frame
  genus_quantity_df <- data.frame(
    Genus = names(genus_quantity_list),
    Quantity = unlist(genus_quantity_list),
    stringsAsFactors = FALSE
  )
  
  return(genus_quantity_df)
}

# Call the function for each stage
genus_quantity_df_stageI <- calculate_genus_quantity(stageI_data, unique_genus_names_stageI, total_samples_stageI)
genus_quantity_df_stageII <- calculate_genus_quantity(stageII_data, unique_genus_names_stageII, total_samples_stageII)
genus_quantity_df_stageIII <- calculate_genus_quantity(stageIII_data, unique_genus_names_stageIII, total_samples_stageIII)
genus_quantity_df_stageIV <- calculate_genus_quantity(stageIV_data, unique_genus_names_stageIV, total_samples_stageIV)



  # Filter out incomplete cases if necessary
genus_quantity_df_stageI <- genus_quantity_df_stageI[complete.cases(genus_quantity_df_stageI), ]
genus_quantity_df_stageII <- genus_quantity_df_stageII[complete.cases(genus_quantity_df_stageII), ]
genus_quantity_df_stageIII <- genus_quantity_df_stageIII[complete.cases(genus_quantity_df_stageIII), ]
genus_quantity_df_stageIV <- genus_quantity_df_stageIV[complete.cases(genus_quantity_df_stageIV), ]



library(ggplot2)

genus_quantity_combined <- rbind(transform(genus_quantity_df_stageI, Stage = "Stage I"),
                                 transform(genus_quantity_df_stageII, Stage = "Stage II"),
                                 transform(genus_quantity_df_stageIII, Stage = "Stage III"),
                                 transform(genus_quantity_df_stageIV, Stage = "Stage IV"))




# Calculate variance across all stages
variance_by_genus <- genus_quantity_combined %>%
  group_by(Genus) %>%
  summarise(variance = var(Quantity, na.rm = TRUE))

# Sort variances in descending order
sorted_variances <- variance_by_genus %>%
  arrange(desc(variance))




filter_sort_and_plot_top_genus <- function(genus_quantity_df, sorted_variances, top_n, filename) {
  # Get the top genus based on variance
  top_n_variances_genus <- head(sorted_variances, top_n)
  top_genus <- top_n_variances_genus$Genus
  
  # Filter data frames for each stage
  filtered_df_stageI <- genus_quantity_df %>%
    filter(Stage == "Stage I" & Genus %in% top_genus)
  
  filtered_df_stageII <- genus_quantity_df %>%
    filter(Stage == "Stage II" & Genus %in% top_genus)
  
  filtered_df_stageIII <- genus_quantity_df %>%
    filter(Stage == "Stage III" & Genus %in% top_genus)
  
  filtered_df_stageIV <- genus_quantity_df %>%
    filter(Stage == "Stage IV" & Genus %in% top_genus)
  
  
  # Combine filtered data frames
  combined_df <- bind_rows(
    mutate(filtered_df_stageI, Stage = "Stage I"),
    mutate(filtered_df_stageII, Stage = "Stage II"),
    mutate(filtered_df_stageIII, Stage = "Stage III"),
    mutate(filtered_df_stageIV, Stage = "Stage IV")
    

  )
  
  
  
  # Plot
  genus_abundance_plot <- ggplot(combined_df, aes(x = factor(Genus, levels = top_genus), y = Quantity, fill = Stage)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    labs(
         x = "Genus", y = "Mean Abundance") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          ) +
    scale_fill_manual(values = c("Stage I" = "#E76BF3", "Stage II" = "#00B0F6", "Stage III" = "#00BF7D", "Stage IV" = "#A3A500"))
  # Save plot to PDF
  pdf(filename)
  print(genus_abundance_plot)
  dev.off()  # Close the PDF device
}



filter_sort_and_plot_top_genus(genus_quantity_combined, sorted_variances, 10, "genus_top10var_final.pdf")















#PLOTTING IN TOP 50 VARIANCES
top_50_variances_genus <- head(sorted_variances, 50)

genus_variances_50 <- ggplot(top_50_variances_genus, aes(x = reorder(Genus, -variance), y = variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Top 50 Genus' with Highest Variances",
       x = "Genus",
       y = "Variance") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())
pdf("genusvarianceplot.pdf")
print(genus_variances_50)
dev.off()  # Close the PDF device
# Plot the values for each stage






#PLOT ON TOTAL GENUS ABUNDANCES - HARD TO READ BECAUSE SO MANY LABELS
combined_plot <- ggplot(genus_quantity_combined, aes(x = Genus, y = Quantity, fill = Stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(title = "", x = "Genus", y = "Total Abundance/Stage Size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Stage I" = "blue", "Stage II" = "red", "Stage III" = "green", "Stage IV" = "purple"))
combined_plot <- combined_plot + theme(axis.text.x = element_text(size = 5))
combined_plot <- combined_plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
combined_plot <- combined_plot + scale_x_discrete(labels = function(x) ifelse(seq_along(x) %% 10 == 0, x, ""))
# Save the combined plot to a PDF file
pdf("genus_abundances_final.pdf")
print(combined_plot)
dev.off()  # Close the PDF device





