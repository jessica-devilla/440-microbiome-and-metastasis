# Print a starting message
cat("Running total abundance analysis by stage...\n")
rm(list = ls(all.names = TRUE))

# load the libraries
suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(ggpubr)
})

source("code/clean_kraken_data.R")

####  load data
#kraken_merge <- readRDS('data/kraken_merge.RDS')
#kraken_merge <- clean_kraken_data(kraken_merge)

kraken_voom_snm <- readRDS('data/kraken_COAD.RDS')
kraken_meta_voom <- readRDS("data/kraken_metaCOAD.RDS")
result_voom <- clean_kraken_data(kraken_voom_snm,kraken_meta_voom)
kraken_data <- result_voom$kraken_data
kraken_meta <- result_voom$kraken_meta

kraken_data$id <- rownames(kraken_data)
kraken_meta$id <- rownames(kraken_meta)


kraken_merge <- merge(kraken_data,kraken_meta[,c('id','pathologic_stage_label','sample_type')], by='id', all.x = TRUE)


stages <- unique(kraken_merge$pathologic_stage_label)

# initialize empty dataframe
total_abund <- data.frame(pathologic_stage_label = stages,
                          mean_val = numeric(length(stages)),
                          std = numeric(length(stages)),
                          n = numeric(length(stages)))

# initialize empty dataframe
indiv_abund <- data.frame(pathologic_stage_label = character(),
                          id = character(),
                          indiv_sums = numeric())

for (i in seq_along(stages)) {
  
  stage_df <- subset(kraken_merge, pathologic_stage_label == stages[i])
  
  sum_cols <- stage_df[, !names(stage_df) %in% c("id", 'pathologic_stage_label', "sample_type")]
  
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
    pathologic_stage_label = rep(stages[i], length(sum_values)),
    id = rownames(stage_df), # Assuming 'id' or rownames can serve as individual IDs
    indiv_sums = sum_values
  )
  
  indiv_abund<- rbind(indiv_abund, temp_df)
}

# Calculate upper and lower limits for error bars
total_abund$upper <- total_abund$mean_val + total_abund$std
total_abund$lower <- total_abund$mean_val - total_abund$std


# Plot the sum of values for each stage category as a bar graph
p <- ggplot(total_abund, aes(x = pathologic_stage_label, y = mean_val, fill = pathologic_stage_label)) +
  geom_bar(stat = "identity") +
  geom_errorbar(data = total_abund, aes(x = pathologic_stage_label, ymin = lower, ymax = upper), width = 0.2, linewidth=0.75,color = "black") +
  labs(x = "Stage", y = "Total Bacterial Abundance") +
  theme_minimal() + ylim(0,5000)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("blue", "#8070FE", "#EAB606","#FC4703"))+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12)) 


print(p)
ggsave(path = "figures", filename = "total_bacterial_abundance_hist.png", bg='white')


p <- ggplot() +
  geom_jitter(data = indiv_abund, aes(x = pathologic_stage_label, y = indiv_sums, color = pathologic_stage_label), width = 0.2, size = 2, alpha = 0.6) +
  geom_point(data = total_abund, aes(x = pathologic_stage_label, y = mean_val), stat = "summary", fun = "mean", size = 3, shape = 19, color = "black") + # Plot mean
  geom_errorbar(data = total_abund, aes(x = pathologic_stage_label, ymin = lower, ymax = upper), width = 0.2, linewidth=0.75,color = "black") +
  labs(x = "Stage", y = "Total Bacterial Abundance") +
  theme_minimal() +
  theme(legend.position = "none")+
  scale_color_manual(values = c("blue", "#8070FE", "#EAB606","#FC4703"))+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12)) 

print(p)
ggsave(path = "figures", filename = "total_bacterial_abundance_jitter.png", bg='white')

ld <- layer_data(last_plot())
(head(ld))

# Calculate significance
# Create a vector to store p-values
p_values <- matrix(NA, nrow = length(unique(indiv_abund$pathologic_stage_label)), ncol = length(unique(indiv_abund$pathologic_stage_label)))

# Perform Wilcoxon rank-sum test for all pairs of stages

for (i in 1:(length(stages) - 1)) {
  for (j in (i + 1):length(stages)) {
    group1 <- indiv_abund$indiv_sums[indiv_abund$pathologic_stage_label == stages[i]]
    group2 <- indiv_abund$indiv_sums[indiv_abund$pathologic_stage_label == stages[j]]
    
    p_values[i, j] <- wilcox.test(group1, group2)$p.value
  }
}

# Fill the lower triangular part of the matrix with the mirrored p-values
p_values[lower.tri(p_values)] <- t(p_values)[lower.tri(p_values)]

# Print or view the p-values matrix
print(p_values)
