#track 
kraken_metaCOAD_subset <- subset(kraken_metaCOAD, sample_size == "primary tumor", select = c("...1", "pathologic_stage_label"))
merged_data_COAD <- merge(kraken_metaCOAD_subset, kraken_COAD, by = "...1")


standardize_stage <- function(stage) {
  # Convert the stage names to lowercase for consistency
  stage <- tolower(stage)
  
  # Map the stage names to their simplified versions
  stage <- ifelse(grepl("stage iv", stage), "Stage IV",
                  ifelse(grepl("stage iii", stage), "Stage III",
                         ifelse(grepl("stage ii", stage), "Stage II",
                                ifelse(grepl("stage i", stage), "Stage I", stage))))
  
  return(stage)
}

merged_data_COAD$pathologic_stage_label <- standardize_stage(merged_data_COAD$pathologic_stage_label)
merged_data_COAD <- merged_data_COAD %>%
  filter(pathologic_stage_label != "not available")


columns_to_modify <- grep("\\.g__", colnames(merged_data_COAD))
colnames(merged_data_COAD)[columns_to_modify] <- sub(".*\\.g__", "", colnames(merged_data_COAD)[columns_to_modify])

library(ggplot2)

create_genus_boxplot <- function(data, genus_column, genus_name, filename) {
  # Subset the data for the specific genus
  genus_subset <- data[, c("...1", "pathologic_stage_label", genus_column)]
  
  # Create the box plot
  box_plot <- ggplot(genus_subset, aes(x = pathologic_stage_label, y = !!rlang::sym(genus_column))) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +  # Add mean points
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "blue") +  # Add standard deviation error bars
    labs(title = paste("Abundance of", genus_name, "Across Stages"),
         x = "Stage",
         y = paste("Abundance of", genus_name)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save plot to PDF
  pdf(filename)
  print(box_plot)
  dev.off()  # Close the PDF device
}

# Usage example:
create_genus_boxplot(merged_data_COAD, "Faecalibacterium", "Faecalibacterium", "faecalibacterium_f_boxplot.pdf")
create_genus_boxplot(merged_data_COAD, "Leptotrichia", "Leptotrichia", "leptotrichia_f_boxplot.pdf")
create_genus_boxplot(merged_data_COAD, "Dorea", "Dorea", "Dorea_f_boxplot.pdf")
create_genus_boxplot(merged_data_COAD, "Odoribacter", "Odoribacter", "Odoribacter_f_boxplot.pdf")
create_genus_boxplot(merged_data_COAD, "Lachnoanaerobaculum", "Lachnoanaerobaculum", "Lachnoanaerobaculum_f_boxplot.pdf")



plot_faecalibacterium_abundance <- function(data, filename) {
  # Aggregate data to calculate total abundance of Faecalibacterium in each stage
  faecalibacterium_aggregated <- data %>%
    group_by(pathologic_stage_label) %>%
    summarise(total_abundance = sum(Faecalibacterium), sample_count = n_distinct(...1)) %>%
    mutate(abundance_per_sample = total_abundance / sample_count)
  
  # Create bar plot
  bar_plot <- ggplot(faecalibacterium_aggregated, aes(x = pathologic_stage_label, y = abundance_per_sample)) +
    geom_bar(stat = "identity", fill = "blue") +
    labs(title = "Total Abundance of Faecalibacterium Divided by Sample Count in Each Stage",
         x = "Stage",
         y = "Abundance per Sample") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save plot to PDF
  pdf(filename)
  print(bar_plot)
  dev.off()  # Close the PDF device
}
plot_faecalibacterium_abundance(merged_data_COAD, "faecalibacterium_abundance_per_sample.pdf")





#cleaningdata
subset_and_remove_column <- function(df, stage_labels) {
  subset_data <- subset(df, pathologic_stage_label %in% stage_labels)
  subset_data <- subset_data[, !names(subset_data) %in% "pathologic_stage_label"]
  return(subset_data)
}
stageI_COADdata <- subset_and_remove_column(merged_data_COAD, c("Stage IA", "Stage IB", "Stage I"))
stageII_COADdata <- subset_and_remove_column(merged_data_COAD, c("Stage IIA", "Stage IIB", "Stage II"))
stageIII_COADdata <- subset_and_remove_column(merged_data_COAD, c("Stage IIIA", "Stage IIIB", "Stage III"))
stageIV_COADdata <- subset_and_remove_column(merged_data_COAD, c("Stage IVA", "Stage IVB", "Stage IV"))



