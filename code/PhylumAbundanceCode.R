install.packages("dplyr")
library("dplyr")
colon_metadata <- Metadata_TCGA_Kraken_17625_Samples %>%
  filter(primary_site == "Colorectal")

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

phylum_names_stageI <- gsub("^.*p__([^\\.]+)\\..*", "\\1", names(stageI_data)[-1])
phylum_names_stageI <- phylum_names_stageI[!grepl("^k__", phylum_names_stageI)]
phylum_names_stageI <- phylum_names_stageI[!grepl("^contaminant", phylum_names_stageI)]

phylum_names_stageII <- gsub("^.*p__([^\\.]+)\\..*", "\\1", names(stageII_data)[-1])
phylum_names_stageII <- phylum_names_stageII[!grepl("^k__", phylum_names_stageII)]
phylum_names_stageII <- phylum_names_stageII[!grepl("^contaminant", phylum_names_stageII)]


phylum_names_stageIII <- gsub("^.*p__([^\\.]+)\\..*", "\\1", names(stageIII_data)[-1])
phylum_names_stageIII <- phylum_names_stageI[!grepl("^k__", phylum_names_stageIII)]
phylum_names_stageIII <- phylum_names_stageI[!grepl("^contaminant", phylum_names_stageIII)]

phylum_names_stageIV <- gsub("^.*p__([^\\.]+)\\..*", "\\1", names(stageIV_data)[-1])
phylum_names_stageIV <- phylum_names_stageI[!grepl("^k__", phylum_names_stageIV)]
phylum_names_stageIV <- phylum_names_stageI[!grepl("^contaminant", phylum_names_stageIV)]

total_samples_stageI <- ncol(stageI_data) - 1  # Subtract 1 to exclude the sample ID column

# Initialize an empty list to store quantities for each phylum
phylum_quantity_list_stageI <- vector("list", length = length(phylum_names_stageI))

# Iterate over each unique phylum name
for (i in seq_along(phylum_names_stageI)) {
  # Filter columns corresponding to the current phylum
  phylum_columns_stageI <- grep(paste0("\\.p__", phylum_names_stageI[i], "\\."), names(stageI_data)[-1], value = TRUE, ignore.case = TRUE)
  
  # Calculate the total amount for the current phylum across all samples
  total_amount_phylum_stageI <- sum(stageI_data[, phylum_columns_stageI])
  
  # Normalize the total amount by the total number of samples in Stage I
  normalized_quantity_stageI <- total_amount_phylum_stageI / total_samples_stageI
  
  # Store the normalized quantity in the list under the phylum name
  phylum_quantity_list_stageI[[phylum_names_stageI[i]]] <- normalized_quantity_stageI
}

# Combine quantities for multiple occurrences of the same phylum
combined_quantities_stageI <- tapply(unlist(phylum_quantity_list_stageI), names(unlist(phylum_quantity_list_stageI)), sum)

# Create a dataframe with phylum names and combined total quantities
phylum_quantity_stageI <- data.frame(Phylum = names(combined_quantities_stageI), Quantity = combined_quantities_stageI)






total_samples_stageII <- ncol(stageII_data) - 1  # Subtract 1 to exclude the sample ID column

# Initialize an empty list to store quantities for each phylum
phylum_quantity_list_stageII <- vector("list", length = length(phylum_names_stageII))

# Iterate over each unique phylum name
for (i in seq_along(phylum_names_stageII)) {
  # Filter columns corresponding to the current phylum
  phylum_columns_stageII <- grep(paste0("\\.p__", phylum_names_stageII[i], "\\."), names(stageII_data)[-1], value = TRUE, ignore.case = TRUE)
  
  # Calculate the total amount for the current phylum across all samples
  total_amount_phylum_stageII <- sum(stageII_data[, phylum_columns_stageII])
  
  # Normalize the total amount by the total number of samples in Stage II
  normalized_quantity_stageII <- total_amount_phylum_stageII / total_samples_stageII
  
  # Store the normalized quantity in the list under the phylum name
  phylum_quantity_list_stageII[[phylum_names_stageII[i]]] <- normalized_quantity_stageII
}

# Combine quantities for multiple occurrences of the same phylum
combined_quantities_stageII <- tapply(unlist(phylum_quantity_list_stageII), names(unlist(phylum_quantity_list_stageII)), sum)

# Create a dataframe with phylum names and combined total quantities
phylum_quantity_stageII <- data.frame(Phylum = names(combined_quantities_stageII), Quantity = combined_quantities_stageII)









total_samples_stageIII <- ncol(stageIII_data) - 1  # Subtract 1 to exclude the sample ID column

# Initialize an empty list to store quantities for each phylum
phylum_quantity_list_stageIII <- vector("list", length = length(phylum_names_stageIII))

# Iterate over each unique phylum name
for (i in seq_along(phylum_names_stageIII)) {
  # Filter columns corresponding to the current phylum
  phylum_columns_stageIII <- grep(paste0("\\.p__", phylum_names_stageIII[i], "\\."), names(stageIII_data)[-1], value = TRUE, ignore.case = TRUE)
  
  # Calculate the total amount for the current phylum across all samples
  total_amount_phylum_stageIII <- sum(stageIII_data[, phylum_columns_stageIII])
  
  # Normalize the total amount by the total number of samples in Stage III
  normalized_quantity_stageIII <- total_amount_phylum_stageIII / total_samples_stageIII
  
  # Store the normalized quantity in the list under the phylum name
  phylum_quantity_list_stageIII[[phylum_names_stageIII[i]]] <- normalized_quantity_stageIII
}

# Combine quantities for multiple occurrences of the same phylum
combined_quantities_stageIII <- tapply(unlist(phylum_quantity_list_stageIII), names(unlist(phylum_quantity_list_stageIII)), sum)

# Create a dataframe with phylum names and combined total quantities
phylum_quantity_stageIII <- data.frame(Phylum = names(combined_quantities_stageIII), Quantity = combined_quantities_stageIII)












total_samples_stageIV <- ncol(stageIV_data) - 1  # Subtract 1 to exclude the sample ID column

# Initialize an empty list to store quantities for each phylum
phylum_quantity_list_stageIV <- vector("list", length = length(phylum_names_stageIV))

# Iterate over each unique phylum name
for (i in seq_along(phylum_names_stageIV)) {
  # Filter columns corresponding to the current phylum
  phylum_columns_stageIV <- grep(paste0("\\.p__", phylum_names_stageIV[i], "\\."), names(stageIV_data)[-1], value = TRUE, ignore.case = TRUE)
  
  # Calculate the total amount for the current phylum across all samples
  total_amount_phylum_stageIV <- sum(stageIV_data[, phylum_columns_stageIV])
  
  # Normalize the total amount by the total number of samples in Stage IV
  normalized_quantity_stageIV <- total_amount_phylum_stageIV / total_samples_stageIV
  
  # Store the normalized quantity in the list under the phylum name
  phylum_quantity_list_stageIV[[phylum_names_stageIV[i]]] <- normalized_quantity_stageIV
}

# Combine quantities for multiple occurrences of the same phylum
combined_quantities_stageIV <- tapply(unlist(phylum_quantity_list_stageIV), names(unlist(phylum_quantity_list_stageIV)), sum)

# Create a dataframe with phylum names and combined total quantities
phylum_quantity_stageIV <- data.frame(Phylum = names(combined_quantities_stageIV), Quantity = combined_quantities_stageIV)

library(ggplot2)

# Combine all data frames into a single data frame with a 'Stage' column
phylum_quantity_combined <- rbind(transform(phylum_quantity_stageI, Stage = "Stage I"),
                                  transform(phylum_quantity_stageII, Stage = "Stage II"),
                                  transform(phylum_quantity_stageIII, Stage = "Stage III"),
                                  transform(phylum_quantity_stageIV, Stage = "Stage IV"))

# Plot combined data without facets
combined_plot <- ggplot(phylum_quantity_combined, aes(x = Phylum, y = Quantity, fill = Stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  labs(title = "Phylum Abundances in the Tumor Microbiome of Colorectal Cancer Samples Across Stages", x = "Phylum", y = "Total Abundance/Stage Size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Stage I" = "blue", "Stage II" = "red", "Stage III" = "green", "Stage IV" = "purple"))

# Save the combined plot to a PDF file
pdf("phylum_abundances_final.pdf")
print(combined_plot)
dev.off()  # Close the PDF device
