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

phylum_names_stageI <- gsub("^.*p__([^\\.]+)\\..*", "\\1", names(stageI_data)[-1])
phylum_names_stageI <- phylum_names_stageI[!grepl("^k__", phylum_names_stageI)]
phylum_names_stageI <- phylum_names_stageI[!grepl("^contaminant", phylum_names_stageI)]
unique_phylum_names_stageI <- unique(phylum_names_stageI)



phylum_names_stageII <- gsub("^.*p__([^\\.]+)\\..*", "\\1", names(stageII_data)[-1])
phylum_names_stageII <- phylum_names_stageII[!grepl("^k__", phylum_names_stageII)]
phylum_names_stageII <- phylum_names_stageII[!grepl("^contaminant", phylum_names_stageII)]
unique_phylum_names_stageII <- unique(phylum_names_stageII)





phylum_names_stageIII <- gsub("^.*p__([^\\.]+)\\..*", "\\1", names(stageIII_data)[-1])
phylum_names_stageIII <- phylum_names_stageI[!grepl("^k__", phylum_names_stageIII)]
phylum_names_stageIII <- phylum_names_stageI[!grepl("^contaminant", phylum_names_stageIII)]
unique_phylum_names_stageIII <- unique(phylum_names_stageIII)




phylum_names_stageIV <- gsub("^.*p__([^\\.]+)\\..*", "\\1", names(stageIV_data)[-1])
phylum_names_stageIV <- phylum_names_stageI[!grepl("^k__", phylum_names_stageIV)]
phylum_names_stageIV <- phylum_names_stageI[!grepl("^contaminant", phylum_names_stageIV)]
unique_phylum_names_stageIV <- unique(phylum_names_stageIV)


total_samples_stageI <- nrow(stageI_data)

phylum_quantity_df_stageI <- data.frame(Phylum = character(), Quantity = numeric(), stringsAsFactors = FALSE)


# Iterate over each unique phylum name
for (phylum_name in unique_phylum_names_stageI){
  # Filter columns corresponding to the current phylum
  phylum_columns_stageI <- grep(paste0("\\.p__",phylum_name,"\\."), names(stageI_data)[-1], value = TRUE, ignore.case = TRUE)
  
  # Calculate the total amount for the current phylum across all samples
  total_amount_phylum_stageI <- sum(stageI_data[, phylum_columns_stageI])
  
  # Normalize the total amount by the total number of samples in Stage I
  normalized_quantity_stageI <- total_amount_phylum_stageI / total_samples_stageI
  phylum_data <- data.frame(Phylum = phylum_name, Quantity = normalized_quantity_stageI, stringsAsFactors = FALSE)
  

  # Store the normalized quantity in the list under the phylum name
  phylum_quantity_df_stageI <- rbind(phylum_quantity_df_stageI, phylum_data)
  
}




total_samples_stageII <- nrow(stageII_data)

phylum_quantity_df_stageII <- data.frame(Phylum = character(), Quantity = numeric(), stringsAsFactors = FALSE)


# Iterate over each unique phylum name
for (phylum_name in unique_phylum_names_stageII){
  # Filter columns corresponding to the current phylum
  phylum_columns_stageII <- grep(paste0("\\.p__",phylum_name,"\\."), names(stageII_data)[-1], value = TRUE, ignore.case = TRUE)
  
  # Calculate the total amount for the current phylum across all samples
  total_amount_phylum_stageII <- sum(stageII_data[, phylum_columns_stageII])
  
  # Normalize the total amount by the total number of samples in Stage I
  normalized_quantity_stageII <- total_amount_phylum_stageII / total_samples_stageII
  phylum_data <- data.frame(Phylum = phylum_name, Quantity = normalized_quantity_stageII, stringsAsFactors = FALSE)
  
  
  # Store the normalized quantity in the list under the phylum name
  phylum_quantity_df_stageII <- rbind(phylum_quantity_df_stageII, phylum_data)
  
}


total_samples_stageIII <- nrow(stageIII_data)

phylum_quantity_df_stageIII <- data.frame(Phylum = character(), Quantity = numeric(), stringsAsFactors = FALSE)


# Iterate over each unique phylum name
for (phylum_name in unique_phylum_names_stageIII){
  # Filter columns corresponding to the current phylum
  phylum_columns_stageIII <- grep(paste0("\\.p__",phylum_name,"\\."), names(stageIII_data)[-1], value = TRUE, ignore.case = TRUE)
  
  # Calculate the total amount for the current phylum across all samples
  total_amount_phylum_stageIII <- sum(stageIII_data[, phylum_columns_stageIII])
  
  # Normalize the total amount by the total number of samples in Stage I
  normalized_quantity_stageIII <- total_amount_phylum_stageIII / total_samples_stageIII
  phylum_data <- data.frame(Phylum = phylum_name, Quantity = normalized_quantity_stageIII, stringsAsFactors = FALSE)
  
  
  # Store the normalized quantity in the list under the phylum name
  phylum_quantity_df_stageIII <- rbind(phylum_quantity_df_stageIII, phylum_data)
  
}





total_samples_stageIV <- nrow(stageIV_data)

phylum_quantity_df_stageIV <- data.frame(Phylum = character(), Quantity = numeric(), stringsAsFactors = FALSE)


# Iterate over each unique phylum name
for (phylum_name in unique_phylum_names_stageIV){
  # Filter columns corresponding to the current phylum
  phylum_columns_stageIV <- grep(paste0("\\.p__",phylum_name,"\\."), names(stageIV_data)[-1], value = TRUE, ignore.case = TRUE)
  
  # Calculate the total amount for the current phylum across all samples
  total_amount_phylum_stageIV <- sum(stageIV_data[, phylum_columns_stageIV])
  
  # Normalize the total amount by the total number of samples in Stage I
  normalized_quantity_stageIV <- total_amount_phylum_stageIV / total_samples_stageIV
  phylum_data <- data.frame(Phylum = phylum_name, Quantity = normalized_quantity_stageIV, stringsAsFactors = FALSE)
  
  
  # Store the normalized quantity in the list under the phylum name
  phylum_quantity_df_stageIV <- rbind(phylum_quantity_df_stageIV, phylum_data)
  
}




library(ggplot2)

# Combine all data frames into a single data frame with a 'Stage' column
phylum_quantity_combined <- rbind(transform(phylum_quantity_df_stageI, Stage = "Stage I"),
                                  transform(phylum_quantity_df_stageII, Stage = "Stage II"),
                                  transform(phylum_quantity_df_stageIII, Stage = "Stage III"),
                                  transform(phylum_quantity_df_stageIV, Stage = "Stage IV"))

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
