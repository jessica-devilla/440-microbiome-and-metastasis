
rm(list = ls(all.names = TRUE))


suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(phyloseq); packageVersion("phyloseq")
  library(ggplot2)
  library(mia)
  library(scater)
  library(umap)
  library(microbiome)
  library(plyr)
  library(vegan)
  library(reshape2)
})

source("code/clean_kraken_data.R")
source("code/run_norm_func.R")
source("code/make_phyloseq_obj.R")
source("code/phyloseq_beta_diversity.R")
source("code/run_zicoseq.R")
#source("code/run_zicoseq_og_plotting.R")



# import and process datasets of interest

#RAW DATA - NORMS NOT WORKING YET
kraken_meta <- readRDS("data/kraken_metaCOAD.RDS")
kraken_data <- readRDS("data/kraken_COAD_raw.RDS")


# UNC DATA
#kraken_meta <- readRDS("data/kraken_meta_norm_filtered.RDS")
#kraken_data<- readRDS("data/kraken_norm_filtered.RDS") # norm'd
#kraken_data <-readRDS("data/kraken_raw_filtered.RDS") # raw

result <- clean_kraken_data(kraken_data, kraken_meta)
kraken_data <- result$kraken_data
kraken_meta <- result$kraken_meta

kraken_data_t <- t(kraken_data)

# Poore et al voom snm
kraken_voom_snm <- readRDS("data/kraken_COAD.RDS")
kraken_meta_voom <- readRDS("data/kraken_metaCOAD.RDS")
result_voom <- clean_kraken_data(kraken_voom_snm,kraken_meta_voom)
kraken_data_voom <- result_voom$kraken_data

## first run the normalization methods on each dataset
# List of normalization methods
norm_methods <- c( "RLE+", "RLE_poscounts", "TSS", "UQ", "CSS",'CLR_poscounts', "logcpm", "CLR+", "MED", "GMPR","DeSEQ")



# Initialize a list to store normalized data frames
normalized_dataframes <- list()

# Loop through each normalization method
for (method in norm_methods) {
  # Try applying normalization function to kraken_data
  tryCatch({
    # Print error message if normalization fails for a method
    cat(paste("Running method: ", method, "\n"))
    # Apply normalization function
    normalized_data <- norm.func(kraken_data_t, method)
    normalized_data <- as.data.frame(normalized_data)
    
    
    # Store the transposed normalized dataframe in the list with method as the key
    normalized_dataframes[[method]] <- normalized_data
  }, error = function(e) {
    # Print error message if normalization fails for a method
    cat(paste("Error occurred for method", method, ":", conditionMessage(e), "\n"))
  })
}

# add raw data to dataframe
normalized_dataframes[["Raw"]]<-kraken_data

#add poore et al voom snm data to dataframe
normalized_dataframes[["Voom-SNM"]] <- kraken_data_voom

## then iterate through each normalized dataset
# create phyloseq object using data and metadata

# Initialize a list to store phyloseq objects
phyloseq_objects <- list()
distance_matrices <- list()

norm_methods <- c("Raw", "Voom-SNM", "RLE+", "RLE_poscounts", "TSS", "UQ", "CSS",'CLR_poscounts', "logcpm", "CLR+", "MED", "GMPR","DeSEQ")

norm_methods <- c("Voom-SNM")
#for (method in names(normalized_dataframes)

# Loop through each normalized data frame
for (method in norm_methods) {
  cat(paste("Running method: ", method, "\n"))
  # Create phyloseq object using create_phyloseq function
  physeq <- make_phyloseq_object(normalized_dataframes[[method]], kraken_meta)
  # Store the phyloseq object in the list with method as the key
  phyloseq_objects[[method]] <- physeq
  
  filename <- paste0("allsamples_", method)
  
  ### run beta diversity function and make plot
  
  #distance_matrix <- physeq_beta_diversity(physeq, dist_methods = c("bray"), name = filename)
  # save distance matrices in list
  #distance_matrices[[method]] <- distance_matrix$bray
  
  ### run zico seq function and plot
  
  kraken_data <- normalized_dataframes[[method]]
  result <- run_zicoseq(kraken_data, kraken_meta, filename)
  
}

zicoObj <- result$zicoObj
pvals <- zicoObj$p.adj.fdr

# calculate alpha diversity for all normalized and plot

# Calculate mantel statistic via spearman correlation between each distance matrix

norm_methods <- c("Raw", "Voom-SNM", "RLE+", "RLE_poscounts", "logcpm","MED", "TSS", "UQ", "CSS",'CLR_poscounts',  "CLR+", "GMPR")



# Initialize an empty matrix to store Mantel statistics
mantel_matrix <- matrix(NA, nrow = length(norm_methods), ncol = length(norm_methods))
rownames(mantel_matrix) <- norm_methods
colnames(mantel_matrix) <- norm_methods
-=# initialize p value matrix
pval_matrix <- matrix(NA, nrow = length(norm_methods), ncol = length(norm_methods))
rownames(pval_matrix) <- norm_methods
colnames(pval_matrix) <- norm_methods

for (i in 1:length(norm_methods)) {
  for (j in 1:length(norm_methods)) {
    cat(paste("Calculating Mantel correlation of ", norm_methods[i], " by ", norm_methods[j], "\n"))
    mantel_result <- mantel(distance_matrices[[i]], distance_matrices[[j]], method = "spearman", permutations = 3)
    mantel_matrix[i, j] <- mantel_result$statistic
    print(mantel_result$statistic)
    print(mantel_result$signif)
    pval_matrix[i, j] <- mantel_result$signif
  }
}

write.table(pval_matrix, file ="data/mantel_pvals.Rdata")
#pval_matrix <- read.table("data/mantel_pvals.Rdata")

write.table(mantel_matrix, file ="data/mantel_mat.Rdata")
#mantel_matrix <- read.table("data/mantel_mat.Rdata")

# plot lower triangle of matrix as a heatmap

matrix_plt <- mantel_matrix

matrix_plt[upper.tri(matrix_plt)]=NA

mantel_df <- melt(matrix_plt)
#mantel_df <- melt(matrix_plt,id.vars=norm_methods)


# Plot the matrix of Mantel statistics
p <- ggplot(mantel_df, aes(Var1, Var2, fill = value)) +
  geom_tile(color='gray') +
  scale_fill_gradient(low = "white", high = "blue",na.value = "white") +
  labs(title = "Mantel Statistic Matrix", x = "Normalization Method", y = "Normalization Method", fill="Mantel\nStatistic") +
  scale_x_discrete(labels = norm_methods) +
  scale_y_discrete(labels = norm_methods)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),  # Adjust size of x-axis labels
        axis.text.y = element_text(size = 13))
print(p)

ggsave("figures/norm_methods_distance_mantel_corr_heatmap_blue_ordered.pdf", plot = p, width=7, height=7)


#mantel_raw <- mantel(distance_matrices[["Raw"]], distance_matrices[["Raw"]], method = "spearman", permutations = 999)

#mantel_voom_raw <- mantel(distance_matrices[["Voom-SNM"]], distance_matrices[["Raw"]], method = "spearman", permutations = 999)

substring<-"Allosalina"
matching_cols <- grep(substring, colnames(kraken_data), value = TRUE)
print(matching_cols)

## calculate p values using bootstrapping procedure

# Initialize an empty matrix to store p-values for matrix comparison
pval_matrix_diff <- matrix(NA, nrow = length(norm_methods), ncol = length(norm_methods))
rownames(pval_matrix_diff) <- norm_methods
colnames(pval_matrix_diff) <- norm_methods

# Number of permutations
n_permutations_diff <- 999

for (i in 1:length(norm_methods)) {
  for (j in 1:length(norm_methods)) {
    cat(paste("Calculating p-value for comparison between", norm_methods[i], "and", norm_methods[j], "\n"))
    
    # If comparing the same matrices, set p-value to 1
    if (i == j) {
      pval_matrix_diff[i, j] <- 1
      next
    }
    
    # Calculate the difference between the distance matrices
    observed_diff <- abs(as.matrix(distance_matrices[[i]]) - as.matrix(distance_matrices[[j]]))
    
    # Perform permutations and calculate difference for each
    permuted_diffs <- replicate(n_permutations_diff, {
      permuted_distance_matrix <- as.matrix(distance_matrices[[i]])
      permuted_distance_matrix <- permuted_distance_matrix[sample(nrow(permuted_distance_matrix)), sample(ncol(permuted_distance_matrix))]
      
      permuted_diff <- abs(permuted_distance_matrix - as.matrix(distance_matrices[[j]]))
      mean(permuted_diff)
    })
    print(paste("Mean Permuted Diffs:", mean(permuted_diffs)))
    
    # Calculate p-value
    observed_statistic_diff <- mean(observed_diff)
    p_value_diff <- mean(permuted_diffs >= observed_statistic_diff)
    pval_matrix_diff[i, j] <- p_value_diff
    print(paste("Observed:", observed_statistic_diff))
    print(paste("Empirical: ",p_value_diff))
  }
}

write.table(pval_matrix_diff, file ="data/distmat_pvals.Rdata")

