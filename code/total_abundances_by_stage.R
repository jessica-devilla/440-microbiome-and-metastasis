# Print a starting message
cat("Running total abundance analysis by stage...\n")

# load the libraries
suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(ggpubr)
})

####  load data
kraken_merge <- readRDS('data/kraken_merge.RDS')
#kraken_merge <- kraken_merge[ , !(names(kraken_merge) %in% c("id",'sample_type'))]

# initialize empty dataframe
total_abund <- data.frame(stage_category = stages,
                          sum_values = numeric(length(stages)),
                          std = numeric(length(stages)),
                          n = numeric(length(stages)))

for (i in seq_along(stages)) {
  
  stage_df <- subset(kraken_merge, stage_category == stages[i])
  
  sum_cols <- stage_df[, !names(stage_df) %in% c("id", "stage_category", "sample_type")]
  
  # Calculate standard deviation
  std_value <- sd(as.matrix(sum_cols))
  
  # Calculate number of observations
  n_value <- nrow(sum_cols)
  
  sum_values <- sum(as.matrix(sum_cols))/n_value
  
  
  # Store the results in the dataframe
  total_abund$sum_values[i] <- sum_values
  total_abund$std[i] <- std_value
  total_abund$n[i] <- n_value
}


# Plot the sum of values for each stage category as a bar graph
ggplot(total_abund, aes(x = stage_category, y = sum_values)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Stage Category", y = "Total Bacterial Abundance") +
  theme_minimal() + ylim(0,5000)
ggsave(path = "figures", filename = "total_bacterial_abundance_hist.png", bg='white')



