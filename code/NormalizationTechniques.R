#https://github.com/wbb121/Norm-Methods-Comparison/blob/main/helper.R

library(tidyr)

#Filtering/Importing Non-Normalized Data,
kraken_orig_otu <- Kraken_TCGA_Raw_Data_17625_Samples
kraken_orig_otu_df <- as.data.frame(kraken_orig_otu)
kraken_raw_COAD <- kraken_orig_otu_df %>%
  filter(...1 %in%kraken_metaCOAD$...1)
colnames(kraken_raw_COAD) <- sub(".*__(.*)$", "\\1", colnames(kraken_raw_COAD))
kraken_metaCOAD_clean <- kraken_metaCOAD
row.names(kraken_metaCOAD_clean) <- kraken_metaCOAD_clean$...1
kraken_metaCOAD_clean <- subset(kraken_metaCOAD_clean, select = -c(...1))

kraken_raw_clean <- kraken_raw_COAD
row.names(kraken_raw_clean) <- kraken_metaCOAD$...1
kraken_raw_clean <- subset(kraken_raw_clean, select = -c(...1))



install.packages("phyloseq")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Use BiocManager to install ConQuR
BiocManager::install("ConQuR")
BiocManager::install("sva")
BiocManager::install("limma")

library(phyloseq)
library(DESeq2)
library(metagenomeSeq)
library(edgeR)
library(compositions)
library(sva)




batch_prep <- subset(kraken_metaCOAD, select = c("...1", "tissue_source_site_label"))
row.names(batch_prep) <- batch_prep$...1
batch_prep <- subset(batch_prep, select = -c(...1))


norm.func <- function(p1,norm_method){
  # remove the all zero genes
  p1 <- p1[rowSums(p1)>0,]
  # remove samples with only 1 non-zero features
  p1 <- p1[,colSums(p1!=0)>1]
  # make names for genes
  rownames(p1) <- make.names(rownames(p1))
  
  if(norm_method=="DeSEQ"){
    sample_ids <- colnames(p1)
    colData <- data.frame(sampleID = sample_ids)
    dds <- DESeqDataSetFromMatrix(countData = p1,
                                  colData = colData,
                                  design = ~1)
    dds <- estimateSizeFactors(dds)
    normalized_counts <- counts(dds, normalized = TRUE)
    norm_p1 <- as.data.frame(normalized_counts)
    norm_p1 <- as.data.frame(t(norm_p1))
    return(norm_p1)
    
  }
  
  if (norm_method == "RLE+") {
    require(DESeq2)
    p1[p1 == 0] <- 1
    
    # Function for RLE+ normalization
    rle.func <- function(tab) {
      metadata <- data.frame(class = factor(1:ncol(tab)))
      tab_dds <- DESeqDataSetFromMatrix(countData = tab, colData = metadata, design = ~class) 
      tab_rle_dds <- estimateSizeFactors(tab_dds)
      tab_norm <- counts(tab_rle_dds, normalized = TRUE)
      return(as.data.frame(tab_norm))
    }
    
    # RLE+ normalization
    final_p1 <- rle.func(p1)
    final_p1 <- as.data.frame(t(final_p1))
    return(final_p1)
  }
  
  # RLE_poscounts, for samples
  # use the non-zero counts to calculate the geometric mean (type="poscounts")
  # performed normalization on the training data, and then performed addon normalization of the test data onto the training data, 
  # to ensure that the normalization of the training data does not in any sense depend on the testing data.
  if(norm_method=="RLE_poscounts"){
    require(DESeq2)
    # function for RLE with poscounts estimator normalization
    rle.poscounts.func <- function(tab){
      metadata <- data.frame(class=factor(1:ncol(tab)))
      tab_dds <- DESeqDataSetFromMatrix(countData=tab,colData=metadata,design=~class) 
      tab_rle_dds <- estimateSizeFactors(tab_dds,type="poscounts")
      tab_norm <- counts(tab_rle_dds, normalized=TRUE)
      as.data.frame(tab_norm)
    }
    # RLE normalization
    norm_p1 <- rle.poscounts.func(p1)
    # let p2 have the same genes as p1 
    norm_p1 <- as.data.frame(t(norm_p1))
    return(norm_p1)
  }
  
  
  # TSS, for samples
  if(norm_method=="TSS"){
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    norm_p1 <- as.data.frame(t(norm_p1))
    return(norm_p1)
  }
  
  
  if(norm_method=="UQ"){
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/quantile(x[x>0])["75%"]))
    norm_p1 <- as.data.frame(t(norm_p1))
    return(norm_p1) 
    
  }
  
  
  if(norm_method=="CSS"){
    require(metagenomeSeq)
    # function for CSS normalization
    css.func <- function(tab){
      tab_mrexperiment <- newMRexperiment(tab)
      tab_css <- cumNorm(tab_mrexperiment)
      tab_norm <- MRcounts(tab_css, norm=T)
      as.data.frame(tab_norm)
    }
    norm_p1 <- css.func(p1)
    norm_p1 <- as.data.frame(t(norm_p1))
    return(norm_p1) 
    
  }
  
  
  if (norm_method == "TMM") {
    require(edgeR)
    
    # Define function for TMM normalization
    tmm.func <- function(tab) {
      tab_dge <- DGEList(counts = tab)
      tab_tmm_dge <- calcNormFactors(tab_dge)
      tab_norm <- cpm(tab_tmm_dge)
      return(as.matrix(tab_norm))
    }
    
    # Perform TMM normalization
    norm_p1 <- tmm.func(p1)
    norm_p1 <- as.data.frame(t(norm_p1))
    
    # Return the normalized matrix
    return(norm_p1)
  }
  
  if(norm_method=="logcpm"){
    require(edgeR)
    # replaced all the 0 abundances with 1 (need integers as input)
    p1[p1==0] <- 1
    # logcpm normalization
    norm_p1 <- as.data.frame(cpm(p1,log=TRUE))
    norm_p1 <- as.data.frame(t(norm_p1))
    return(norm_p1)
  }
  
  
  if(norm_method == 'rarefy'){
    require(phyloseq)
    ps <- otu_table(as.matrix(p1), taxa_are_rows = FALSE)
    rarefied_ps <- rarefy_even_depth(ps, sample.size = 1)
    norm_p1 <- as.data.frame(otu_table(rarefied_ps))
    norm_p1 <- as.data.frame(t(norm_p1))
    return(norm_p1)
  }
  
  
  if (norm_method == "CLR+") {
    require(compositions)
    
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1, 2, function(x) x / sum(x)))
    
    # Replace all the 0 abundances with 0.65 times minimum non-zero abundance
    min_non_zero <- min(norm_p1[norm_p1 != 0])
    pseudo_count <- min_non_zero * 0.65
    
    norm_p1[norm_p1 == 0] <- pseudo_count
    
    # CLR transformation
    trans_p1 <- as.data.frame(apply(norm_p1, 2, function(x) clr(x)))
    
    # Ensure no negative values
    trans_p1[trans_p1 < 0] <- 0
    trans_p1 <- as.data.frame(t(trans_p1))
    
    return(trans_p1)
  }
  
  
  if (norm_method == "combat") {
    require(sva)
    
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1, 2, function(x) x / sum(x)))
    
    # Log transformation
    norm_p1 <- log(norm_p1 + 1)  # Adding 1 to avoid log(0) which is undefined
    
    # Replace all the 0 abundances with 0.65 times the minimum non-zero abundance
    norm_p1[norm_p1 == 0] <- min(norm_p1[norm_p1 != 0]) * 0.65
    
    # Check for genes with no variation within a single batch
    zero_var_genes <- apply(norm_p1, 1, function(row) all(row == 0))
    norm_p1 <- norm_p1[!zero_var_genes, ]
    
    # Create batch factor (assuming p1 represents your data)
    batch_factor <- as.factor(batch_prep$tissue_source_site_label)
    
    # Apply ComBat normalization
    correct_combat <- ComBat(norm_p1, batch = batch_factor, ref.batch = "Indivumed" )
    correct_combat_df <- as.data.frame(correct_combat)
    
    return(correct_combat_df)
  }
  
  if(norm_method=="MED"){
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/median(x[x>0])))
    # let p2 have the same genes as p1 
    norm_p1 <- as.data.frame(t(norm_p1))
    return(norm_p1)
  }
  
  if (norm_method == "GMPR") {
    require(GUniFrac)
    
    # Function for GMPR normalization
    gmpr.func <- function(tab) {
      gmpr_size_factor <- GMPR(tab)
      tab_norm <- as.data.frame(t(t(tab) / gmpr_size_factor))
      return(tab_norm)
    }
    
    # GMPR normalization
    final_p1 <- gmpr.func(p1)
    final_p1 <- as.data.frame(t(final_p1))
    return(final_p1)
  }
  
  if (norm_method == "CLR_poscounts") {
    require(compositions)
    
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1, 2, function(x) x / sum(x)))
    
    # Calculate the minimum non-zero abundance
    min_non_zero <- min(norm_p1[norm_p1 != 0])
    
    # Add a pseudo-count to the data to avoid zero values
    pseudo_count <- min_non_zero * 0.65
    norm_p1[norm_p1 == 0] <- pseudo_count
    
    # CLR transformation
    trans_p1 <- as.data.frame(apply(norm_p1, 2, function(x) clr(x)))
    
    # Ensure no negative values
    trans_p1[trans_p1 < 0] <- 0
    trans_p1 <- as.data.frame(t(trans_p1))
    
    return(trans_p1)
  }
  
  
}

#Normalization Techniques - Calling Function Above
transposed_df_kraken <- t(kraken_raw_clean)
transposed_df_kraken <- as.data.frame(transposed_df_kraken)
dim(transposed_df_kraken)

tot_TSS <- norm.func(transposed_df_kraken, 'TSS')
tot_MED <- norm.func(transposed_df_kraken, "MED")
tot_GMPR <- norm.func(transposed_df_kraken, 'GMPR')
tot_UQ <- norm.func(transposed_df_kraken, 'UQ')
tot_CSS <- norm.func(transposed_df_kraken, 'CSS')
tot_DEQ <- norm.func(transposed_df_kraken, 'DeSEQ')
tot_RLE <- norm.func(transposed_df_kraken, 'RLE+')
tot_RLEpos <- norm.func(transposed_df_kraken, 'RLE_poscounts')
#tot_TMM <- norm.func(transposed_df_kraken, 'TMM')
tot_logcpm <- norm.func(transposed_df_kraken, 'logcpm')
tot_rarefy <- norm.func(transposed_df_kraken, 'rarefy')
tot_CLR <- norm.func(transposed_df_kraken, 'CLR+')
tot_CLRpos <- norm.func(transposed_df_kraken, 'CLR_poscounts')


row.names(kraken_COAD_genus) <- kraken_COAD_genus$...1
kraken_COAD_genus_clean <- subset(kraken_COAD_genus, select = -c(...1))

data_list <- list(
  TSS = tot_TSS,
  MED = tot_MED,
  GMPR = tot_GMPR,
  UQ = tot_UQ,
  CSS = tot_CSS,
  DEQ = tot_DEQ,
  RLE = tot_RLE,
  RLEpos = tot_RLEpos,
  #TMM = tot_TMM,  # Assuming tot_TMM is defined
  logcpm = tot_logcpm,
  CLR = tot_CLR,
  CLRpos = tot_CLRpos,
  Raw = kraken_raw_clean,
  VoomSNM = kraken_COAD_genus_clean
)

# PREFORM UMAP 
library(umap)
library(purrr)
library(ggplot2)
library(dplyr) 
library(Rtsne)

set.seed(123)
umap_list <- lapply(data_list, function(data) {
  set.seed(123)  # Set seed for reproducibility
  umap(data)
})


library(purrr)

# Create a list of umap results
umap_df_list <- map2(data_list, names(data_list), function(data, method) {
  umap_result <- umap(data)
  data.frame(UMAP1 = umap_result$layout[,1], UMAP2 = umap_result$layout[,2], Method = method)
})

# Combine umap results into a single data frame
umap_df <- bind_rows(umap_df_list)

library(ggplot2)
pdf("umap_plot_total_playing.pdf")
umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Method)) +
  geom_point(size = .5) +
  theme_minimal() +
  labs(title = "UMAP Visualization of Different Normalization Methods")
print(umap_plot)
dev.off()


#PREFORM TSNE
perplexity <- 30  # Set your desired perplexity value
tsne_list <- lapply(data_list, function(data) Rtsne(data, perplexity = perplexity)$Y)
tsne_df <- do.call(rbind, Map(function(data, method) {
  data.frame(V1 = data[,1], V2 = data[,2], Method = method)
}, tsne_list, names(tsne_list)))

# Plot t-SNE visualization
# Plot t-SNE visualization
pdf("tnse_plot_total_playing.pdf", width = 10, height = 10)
tsne_plot <- ggplot(tsne_df, aes(x = V1, y = V2, color = Method)) +
  geom_point(size = 0.4) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15, margin = margin(t = 10), color = "black"),
    axis.title.y = element_text(size = 15, margin = margin(r = 10), color = "black"),
    axis.text.x = element_text(margin = margin(t = 10), size = 10, color = "black"),  # Increase x-axis label size
    axis.text.y = element_text(margin = margin(r = 10), size = 10, color = "black"),
    axis.line = element_line(color = "black"),  # Change color of axis lines
    panel.grid.major = element_blank(),
    legend.text = element_text(size = 10),  # Increase legend text size
    legend.key.size = unit(1, "lines"),  # Increase space between legend elements
    panel.grid.minor = element_blank(), 
    plot.title = element_text(size = 20, margin = margin(b=10), hjust = 0.5),
    legend.title =element_blank()) +
  labs(title = "Normalization Methods T-SNE Visualization")+
  guides(color = guide_legend(key.width = unit(10, "cm"), key.height = unit(10, "cm"))) # Adjust key size here
print(tsne_plot)
dev.off()  # Close the PDF device


data_list_with_stages <- lapply(data_list, function(df) {
  # Match row names with kraken_metaCOAD_clean
  matched_rows <- intersect(rownames(df), rownames(kraken_metaCOAD_clean))

  # Add stage_label column with corresponding stage labels
  df$pathologic_stage_label <- kraken_metaCOAD_clean[matched_rows, "pathologic_stage_label"]  
  # Reorder columns to make stage_label the first column
  df <- df[, c("pathologic_stage_label", setdiff(colnames(df), "pathologic_stage_label"))]
  
  return(df)
})




methods <- unique(names(data_list_with_stages))
print(methods)

kruskal_results <- list()
spearman_results <- list()

for (method in methods) {
  current_df <- data_list_with_stages[[method]]
  current_df$id <- rownames(current_df)  # Add row names as a new column named 'id'
  current_df <- current_df[, c("id", "pathologic_stage_label", setdiff(names(current_df), c("id", "pathologic_stage_label")))]
  
  # Check if the dataframe has the stage_label column
  if ("id" %in% colnames(current_df)) {
    kruskal_name <- paste0("kruskal_", method)
    output_file_kruskal <- paste0(kruskal_name, ".pdf")
    metadata_df <- kraken_metaCOAD[kraken_metaCOAD$id %in% current_df$id, ]
    
    kruskal_result <- perform_kruskal_and_plot_abundance(input_df = current_df,
                                                         metadata_df = metadata_df,
                                                         n_top = 5,
                                                         output_file = output_file_kruskal)
    kruskal_results[[method]] <- kruskal_result
    
    
    assign(kruskal_name, kruskal_result)
    
    spear_name <- paste0("spear_", method)
    output_file_spear <- paste0(spear_name, ".pdf")
    
    spear_result <- perform_spearman_and_plot_abundance(input_df = current_df,
                                                        metadata_df = metadata_df,
                                                        n_top = 5,
                                                        output_file = output_file_spear)
    spearman_results[[method]] <- spear_result
    
    assign(spear_name, spear_result)
    
  } else {
    cat("Dataframe", i, "does not contain the stage_label column.\n")
  }
}

combined_kruskal <- do.call(rbind, kruskal_results)

# Combine Spearman results into a single dataframe
combined_spearman <- do.call(rbind, spearman_results)
combined_spearman$Method <- sub("\\..*", "", rownames(combined_spearman))
combined_spearman <- combined_spearman[!combined_spearman$Method %in% c("RLE", "logcpm"), ]


# Perform rank sum on combined dataframes
ranksum_kruskal_submit <- combined_kruskal %>%
  group_by(taxon) %>%
  summarise(total_rank = sum(rank)) %>%
  arrange(total_rank)

ranksum_spearman_submit <- combined_spearman %>%
  group_by(taxon) %>%
  summarise(total_rank = sum(rank)) %>%
  arrange(total_rank)


top_taxa <- ranksum_spear_submit$taxon[1:10]  # Selecting top 10 taxa

# Filter the combined PCA dataset to include only the top taxa
subset_data <- combined_spearman %>%
  filter(taxon %in% top_taxa)

# Reorder the 'Method' factor variable
subset_data$Method <- factor(subset_data$Method, 
                             levels = c("Raw", "VoomSNM", "DEQ", "GMPR", "UQ", "TSS", 
                                        "RLEpos", "MED", "CLR", "CLRpos", "CSS"))

heatmap_plot <- ggplot(subset_data, aes(x = Method, y = factor(taxon, levels = top_taxa), fill = abs(correlation))) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "#012300", high = "white") +
  theme_minimal() +
  labs(title = "Taxa with Top Absolute Correlation Across Normalization Techniques",
       x = "Normalization Method", y = "Genus", fill = "Correlation") +
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 22, margin = margin(t = 20), color = "black"),
        axis.title.y = element_text(size = 22, margin = margin(r = 15), color = "black"),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 15, hjust = 0.5),  # Increase legend text size
        legend.title = element_text(size = 22, hjust = 0.5),  # Increase legend title size
        legend.position = "right",  # Position
        plot.title = element_text(size = 30, margin = margin(b = 20), hjust = 0.5),
        legend.key.size = unit(2, "lines"),  # Increase space between legend elements
        axis.title = element_text(size = 12))

pdf("heatmap_spear_normtech.pdf", width = 20, height = 10)
print(heatmap_plot)
dev.off()









for (method in names(spearman_results)) {
  current_df <- spearman_results[[method]]
  
  pdf(paste0("spearman_correlation_", method, ".pdf"))  
  top_taxa <- current_df %>%
    slice_head(n = 2)  # Change the number to select the top N taxa
  
  # Convert sign of correlation to a factor with two levels
  current_df$correlation_sign <- ifelse(current_df$correlation >= 0, "positive", "negative")
  
  # Volcano plot with labeled top ranked taxa
  volcano_plot_labeled <- ggplot(current_df, aes(x = correlation, y = -log10(p_adjust), color = correlation_sign)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_text(data = top_taxa, aes(label = taxon), hjust = -0.2, vjust = -0.5, size = 4, color = "black", bg = "transparent") +
    labs(title = bquote(paste("Spearman's ", rho, " between Stage and Genus in ", .(method), " on Tumor Microbiome Data")),
         x = expression(paste("Spearman's ", rho)),
         y = "-log10(Adjusted p-value)",
         color = "Correlation") +
    scale_color_manual(values = c("positive" = "#EAB606", "negative" = "#8070FE"), 
                       labels = c("positive" = "Positive", "negative" = "Negative")) +
    theme(
      plot.title = element_text(size = 10, margin = margin(b=10), hjust = 0.5),
      legend.position = "bottom",
      axis.line = element_line(color = "black")) +  # Add axis lines for both x and y axes
    xlim(c(-0.3, 0.3)) +  # Adjust the limits according to your data range
    theme(panel.background = element_blank())  # Remove plot background
  
  # Print the plot
  print(volcano_plot_labeled)
  
  dev.off()
  
}
  
  
  



