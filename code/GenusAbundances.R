colon_metadata <- kraken_metaCOAD
Kraken_TCGA_Voom_SNM_Plate_Center_Filtering_Data <- kraken_df
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

genus_names_stageI <- gsub("^.*g__([^\\.]+).*", "\\1", names(stageI_data)[-1])
genus_names_stageI <- genus_names_stageI[!grepl("^k__", genus_names_stageI)]
genus_names_stageI <- genus_names_stageI[!grepl("^contaminant", genus_names_stageI)]
unique_genus_names_stageI <- unique(genus_names_stageI)


genus_names_stageII <- gsub("^.*g__([^\\.]+).*", "\\1", names(stageI_data)[-1])
genus_names_stageII <- genus_names_stageII[!grepl("^k__", genus_names_stageII)]
genus_names_stageII <- genus_names_stageII[!grepl("^contaminant", genus_names_stageII)]
unique_genus_names_stageII <- unique(genus_names_stageII)





genus_names_stageIII <- gsub("^.*g__([^\\.]+).*", "\\1", names(stageI_data)[-1])
genus_names_stageIII <- genus_names_stageI[!grepl("^k__", genus_names_stageIII)]
genus_names_stageIII <- genus_names_stageI[!grepl("^contaminant", genus_names_stageIII)]
unique_genus_names_stageIII <- unique(genus_names_stageIII)



genus_names_stageIV <- gsub("^.*g__([^\\.]+).*", "\\1", names(stageI_data)[-1])
genus_names_stageIV <- genus_names_stageI[!grepl("^k__", genus_names_stageIV)]
genus_names_stageIV <- genus_names_stageI[!grepl("^contaminant", genus_names_stageIV)]
unique_genus_names_stageIV <- unique(genus_names_stageIV)


total_samples_stageI <- nrow(stageI_data)

genus_quantity_df_stageI <- data.frame(Genus = character(), Quantity = numeric(), stringsAsFactors = FALSE)


# Iterate over each unique genus name
for (genus_name in unique_genus_names_stageI){
  # Filter columns corresponding to the current genus
  genus_columns_stageI <- grep(paste0("\\.g__",genus_name), names(stageI_data)[-1], value = TRUE, ignore.case = TRUE)
  
  # Calculate the total amount for the current phylum across all samples
  total_amount_genus_stageI <- sum(stageI_data[, genus_columns_stageI])
  
  # Normalize the total amount by the total number of samples in Stage I
  normalized_quantity_stageI <- total_amount_genus_stageI / total_samples_stageI
  genus_data <- data.frame(Genus = genus_name, Quantity = normalized_quantity_stageI, stringsAsFactors = FALSE)
  
  
  # Store the normalized quantity in the list under the phylum name
  genus_quantity_df_stageI <- rbind(genus_quantity_df_stageI, genus_data)
  
}






total_samples_stageII <- nrow(stageII_data)

genus_quantity_df_stageII <- data.frame(Genus = character(), Quantity = numeric(), stringsAsFactors = FALSE)


# Iterate over each unique genus name
for (genus_name in unique_genus_names_stageII){
  # Filter columns corresponding to the current genus
  genus_columns_stageII <- grep(paste0("\\.g__",genus_name), names(stageII_data)[-1], value = TRUE, ignore.case = TRUE)
  
  # Calculate the total amount for the current phylum across all samples
  total_amount_genus_stageII <- sum(stageII_data[, genus_columns_stageII])
  
  # Normalize the total amount by the total number of samples in Stage I
  normalized_quantity_stageII <- total_amount_genus_stageII / total_samples_stageII
  genus_data <- data.frame(Genus = genus_name, Quantity = normalized_quantity_stageII, stringsAsFactors = FALSE)
  
  
  # Store the normalized quantity in the list under the phylum name
  genus_quantity_df_stageII <- rbind(genus_quantity_df_stageII, genus_data)
  
}




total_samples_stageIII <- nrow(stageIII_data)

genus_quantity_df_stageIII <- data.frame(Genus = character(), Quantity = numeric(), stringsAsFactors = FALSE)


# Iterate over each unique genus name
for (genus_name in unique_genus_names_stageIII){
  # Filter columns corresponding to the current genus
  genus_columns_stageIII <- grep(paste0("\\.g__",genus_name), names(stageIII_data)[-1], value = TRUE, ignore.case = TRUE)
  
  # Calculate the total amount for the current phylum across all samples
  total_amount_genus_stageIII <- sum(stageIII_data[, genus_columns_stageIII])
  
  # Normalize the total amount by the total number of samples in Stage I
  normalized_quantity_stageIII <- total_amount_genus_stageIII / total_samples_stageIII
  genus_data <- data.frame(Genus = genus_name, Quantity = normalized_quantity_stageIII, stringsAsFactors = FALSE)
  
  
  # Store the normalized quantity in the list under the phylum name
  genus_quantity_df_stageIII <- rbind(genus_quantity_df_stageIII, genus_data)
  
}










total_samples_stageIV <- nrow(stageIV_data)

genus_quantity_df_stageIV <- data.frame(Genus = character(), Quantity = numeric(), stringsAsFactors = FALSE)


# Iterate over each unique genus name
for (genus_name in unique_genus_names_stageIV){
  # Filter columns corresponding to the current genus
  genus_columns_stageIV <- grep(paste0("\\.g__",genus_name), names(stageIV_data)[-1], value = TRUE, ignore.case = TRUE)
  
  # Calculate the total amount for the current phylum across all samples
  total_amount_genus_stageIV <- sum(stageIV_data[, genus_columns_stageIV])
  
  # Normalize the total amount by the total number of samples in Stage I
  normalized_quantity_stageIV <- total_amount_genus_stageIV / total_samples_stageIV
  genus_data <- data.frame(Genus = genus_name, Quantity = normalized_quantity_stageIV, stringsAsFactors = FALSE)
  
  
  # Store the normalized quantity in the list under the phylum name
  genus_quantity_df_stageIV <- rbind(genus_quantity_df_stageIV, genus_data)
  
}
genus_quantity_df_stageIV <- genus_quantity_df_stageIV[complete.cases(genus_quantity_df_stageIV), ]
genus_quantity_df_stageIII <- genus_quantity_df_stageIII[complete.cases(genus_quantity_df_stageIII), ]
genus_quantity_df_stageII <- genus_quantity_df_stageII[complete.cases(genus_quantity_df_stageII), ]
genus_quantity_df_stageI <- genus_quantity_df_stageI[complete.cases(genus_quantity_df_stageI), ]





#PLOTTING BASED ON VARIANCE

library(ggplot2)

# Combine all data frames into a single data frame with a 'Stage' column
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


#PLOTTING IN TOP 50 VARIANCES
top_50_variances_genus <- head(sorted_variances, 50)

top_50_genus <- top_50_variances_genus$Genus

filtered_df_stageI <- genus_quantity_df_stageI %>%
  filter(Genus %in% top_50_genus)

filtered_df_stageII <- genus_quantity_df_stageII %>%
  filter(Genus %in% top_50_genus)

filtered_df_stageIII <- genus_quantity_df_stageIII %>%
  filter(Genus %in% top_50_genus)

filtered_df_stageIV <- genus_quantity_df_stageIV %>%
  filter(Genus %in% top_50_genus)

# Combine filtered data frames
combined_top50vargenus <- bind_rows(
  mutate(filtered_df_stageI, Stage = "Stage I"),
  mutate(filtered_df_stageII, Stage = "Stage II"),
  mutate(filtered_df_stageIII, Stage = "Stage III"),
  mutate(filtered_df_stageIV, Stage = "Stage IV")
)


genusabudance_top50variance <- ggplot(combined_top50vargenus, aes(x = factor(Genus,levels = top_50_variances_genus$Genus), y = Quantity, fill = Stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(title = "Values of Top 50 Genus' with Highest Variances across Stages",
       x = "Genus", y = "Total Abundance/Stage Size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Stage I" = "blue", "Stage II" = "red", "Stage III" = "green", "Stage IV" = "purple"))


pdf("genus_top50var_final.pdf")
print(genusabudance_top50variance)
dev.off()  # Close the PDF device
# Plot the values for each stage




#PLOTTING TOP 25 VARIANCE
top_25_variances_genus <- head(sorted_variances, 25)

top_25_genus <- top_25_variances_genus$Genus

filteredTOP25_df_stageI <- genus_quantity_df_stageI %>%
  filter(Genus %in% top_25_genus)

filteredTOP25_df_stageII <- genus_quantity_df_stageII %>%
  filter(Genus %in% top_25_genus)

filteredTOP25_df_stageIII <- genus_quantity_df_stageIII %>%
  filter(Genus %in% top_25_genus)

filteredTOP25_df_stageIV <- genus_quantity_df_stageIV %>%
  filter(Genus %in% top_25_genus)

# Combine filtered data frames
combined_top25vargenus <- bind_rows(
  mutate(filteredTOP25_df_stageI, Stage = "Stage I"),
  mutate(filteredTOP25_df_stageII, Stage = "Stage II"),
  mutate(filteredTOP25_df_stageIII, Stage = "Stage III"),
  mutate(filteredTOP25_df_stageIV, Stage = "Stage IV")
)


genusabudance_top25variance <- ggplot(combined_top25vargenus, aes(x = factor(Genus,levels = top_25_variances_genus$Genus), y = Quantity, fill = Stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(title = "Values of Top 25 Genus' with Highest Variances across Stages",
       x = "Genus", y = "Total Abundance/Stage Size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Stage I" = "blue", "Stage II" = "red", "Stage III" = "green", "Stage IV" = "purple"))


pdf("genus_top25var_final.pdf")
print(genusabudance_top25variance)
dev.off()  # Close the PDF device
# Plot the values for each stage




























#PLOT ON TOTAL GENUS ABUNDANCES - HARD TO READ BECAUSE SO MANY LABELS
combined_plot <- ggplot(genus_quantity_combined, aes(x = Genus, y = Quantity, fill = Stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(title = "Average Genus Abundances in a Tumor Microbiome of Colorectal Cancer Across Stages", x = "Genus", y = "Total Abundance/Stage Size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Stage I" = "blue", "Stage II" = "red", "Stage III" = "green", "Stage IV" = "purple"))
combined_plot <- combined_plot + theme(axis.text.x = element_text(size = 5))
combined_plot <- combined_plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
combined_plot <- combined_plot + scale_x_discrete(labels = function(x) ifelse(seq_along(x) %% 10 == 0, x, ""))


# Save the combined plot to a PDF file
pdf("genus_abundances_final.pdf")
print(combined_plot)
dev.off()  # Close the PDF device








#PLOT ON GENUS DIFFERENCE BETWEEN STAGES - PLOTTED TOP 50 
genus_difference <- genus_quantity_combined %>%
  group_by(Genus) %>%
  summarize(Difference = max(Quantity) - min(Quantity)) %>%
  arrange(desc(Difference))

# Select the top 100 genera with the largest difference
top_50_genus_diff <- head(genus_difference$Genus, 50)

# Filter the dataframe to include only the top 100 genera
top_50_genus_diff_data <- genus_quantity_combined %>%
  filter(Genus %in% top_50_genus_diff)

# Plot the data for the top 100 genera
top_50_plot_genus_diff <- ggplot(top_50_genus_diff_data, aes(x = Genus, y = Quantity, fill = Stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(title = "Top 50 Genus Abundances with Largest Differences Between Strains",
       x = "Genus", y = "Total Abundance/Stage Size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Stage I" = "blue", "Stage II" = "red", "Stage III" = "green", "Stage IV" = "purple"))

# Save the plot to a PDF file
pdf("top_50_genus_abundances_differences.pdf")
print(top_50_plot_genus_diff)
dev.off()  # Close the PDF device

#possibly wilcox test#
#better way to differentiate between the changes in genus then just differences#
#maybe look at the change between stage 2 and stage 3, then stage 3 and stage 4#
#track the pattern of certain strains across stages and see if this is consistent across all cells#





