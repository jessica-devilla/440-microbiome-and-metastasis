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

stages <- unique(kraken_merge$stage_category)

# initialize empty dataframe
total_abund <- data.frame(stage_category = stages,
                          mean_val = numeric(length(stages)),
                          std = numeric(length(stages)),
                          n = numeric(length(stages)))

# initialize empty dataframe
indiv_abund <- data.frame(stage_category = character(),
                          id = character(),
                          indiv_sums = numeric())

for (i in seq_along(stages)) {
  
  stage_df <- subset(kraken_merge, stage_category == stages[i])
  
  sum_cols <- stage_df[, !names(stage_df) %in% c("id", "stage_category", "sample_type")]
  
  sum_values <- rowSums(as.matrix(sum_cols))
  
  mean_value <- mean(sum_values)
  
  # Calculate standard deviation
  std_value <- sd(sum_values)
  
  # Calculate number of observations
  n_value <- length(sum_values)
  
  
  # Store the results in the dataframe
  total_abund$mean_val[i] <- mean_value
  total_abund$std[i] <- std_value
  total_abund$n[i] <- n_value
  
  temp_df <- data.frame(
    stage_category = rep(stages[i], length(sum_values)),
    id = rownames(stage_df), # Assuming 'id' or rownames can serve as individual IDs
    indiv_sums = sum_values
  )
  
  indiv_abund<- rbind(indiv_abund, temp_df)
}

# Calculate upper and lower limits for error bars
total_abund$upper <- total_abund$mean_val + total_abund$std
total_abund$lower <- total_abund$mean_val - total_abund$std


# Plot the sum of values for each stage category as a bar graph
p <- ggplot(total_abund, aes(x = stage_category, y = mean_val, fill = stage_category)) +
  geom_bar(stat = "identity") +
  geom_errorbar(data = total_abund, aes(x = stage_category, ymin = lower, ymax = upper), width = 0.2, color = "black") +
  labs(x = "Stage Category", y = "Total Bacterial Abundance") +
  theme_minimal() + ylim(0,5000)+
  theme(legend.position = "none")
print(p)
ggsave(path = "figures", filename = "total_bacterial_abundance_hist.png", bg='white')


p <- ggplot() +
  geom_jitter(data = indiv_abund, aes(x = stage_category, y = indiv_sums, color = stage_category), width = 0.2, size = 2, alpha = 0.6) +
  geom_point(data = total_abund, aes(x = stage_category, y = mean_val), stat = "summary", fun = "mean", size = 3, shape = 19, color = "black") + # Plot mean
  geom_errorbar(data = total_abund, aes(x = stage_category, ymin = lower, ymax = upper), width = 0.2, color = "black") +
  labs(x = "Stage Category", y = "Total Bacterial Abundance") +
  theme_minimal() +
  theme(legend.position = "none")

print(p)
ggsave(path = "figures", filename = "total_bacterial_abundance_jitter.png", bg='white')

ld <- layer_data(last_plot())
(head(ld))

# Calculate significance
# Create a vector to store p-values
p_values <- matrix(NA, nrow = length(unique(indiv_abund$stage_category)), ncol = length(unique(indiv_abund$stage_category)))

# Perform Wilcoxon rank-sum test for all pairs of stages

for (i in 1:(length(stages) - 1)) {
  for (j in (i + 1):length(stages)) {
    group1 <- indiv_abund$indiv_sums[indiv_abund$stage_category == stages[i]]
    group2 <- indiv_abund$indiv_sums[indiv_abund$stage_category == stages[j]]
    
    p_values[i, j] <- wilcox.test(group1, group2)$p.value
  }
}

# Fill the lower triangular part of the matrix with the mirrored p-values
p_values[lower.tri(p_values)] <- t(p_values)[lower.tri(p_values)]

# Print or view the p-values matrix
print(p_values)
