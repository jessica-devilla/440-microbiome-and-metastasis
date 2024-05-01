#https://github.com/wbb121/Norm-Methods-Comparison/blob/main/helper.R



install.packages("phyloseq")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Use BiocManager to install ConQuR
BiocManager::install("ConQuR")
BiocManager::install("sva")
BiocManager::install("limma")
BiocManager::install("DESeq2")

library(phyloseq)
library(DESeq2)



batch_prep <- subset(kraken_metaCOAD_RNASeq_IlluminaGA_UNC, select = c("...1", "tissue_source_site_label"))
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
    final_p1 <- norm_p1
    return(final_p1)
  }
  
  
  # TSS, for samples
  if(norm_method=="TSS"){
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    return(norm_p1)
  }
  
  
  if(norm_method=="UQ"){
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/quantile(x[x>0])["75%"]))
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
    
    # Return the normalized matrix
    return(norm_p1)
  }
  
  if(norm_method=="logcpm"){
    require(edgeR)
    # replaced all the 0 abundances with 1 (need integers as input)
    p1[p1==0] <- 1
    # logcpm normalization
    norm_p1 <- as.data.frame(cpm(p1,log=TRUE))
    return(norm_p1)
  }
  
  
  if(norm_method == 'rarefy'){
    require(phyloseq)
    ps <- otu_table(as.matrix(p1), taxa_are_rows = FALSE)
    rarefied_ps <- rarefy_even_depth(ps, sample.size = 1)
    norm_p1 <- as.data.frame(otu_table(rarefied_ps))
    return(norm_p1)
  }
  
  
  # CLR+, for samples
  # based on TSS normalized data, add a pseudo count 0.65*minimum to zero values
  if(norm_method == "CLR+") {
    require(compositions)
    norm_p1 <- as.data.frame(apply(p1, 2, function(x) x / sum(x)))
    
    # Replace all the 0 abundances with 0.65 times minimum non-zero abundance
    epsilon <- 1e-6  # small constant to avoid taking log of zero
    min_non_zero <- min(norm_p1[norm_p1 != 0])
    norm_p1[norm_p1 == 0] <- min_non_zero * 0.65
    trans_p1 <- as.data.frame(apply(norm_p1, 2, function(x) clr(x + epsilon)))
    offset <- 3  
    trans_p1 <- trans_p1 + offset
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
    return(final_p1)
  }
  
  if(norm_method=="CLR_poscounts"){
    require(compositions)
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    # clr transformation
    trans_p1 <- as.data.frame(apply(norm_p1,2,function(x) clr(x)))
    # let p2 have the same genes as p1 
    return (trans_p1)
  }
  
  
}

#Normalization Techniques - Calling Function Above
transposed_df_krakenILUNC <- t(kraken_raw_COADRNA_IlluminaGA_UNC_clean)

ILUNC_TSS <- norm.func(transposed_df_krakenILUNC, 'TSS')
ILUNC_MED <- norm.func(transposed_df_krakenILUNC, 'MED')
ILUNC_GMPR <- norm.func(transposed_df_krakenILUNC, 'GMPR')
ILUNC_UQ <- norm.func(transposed_df_krakenILUNC, 'UQ')
ILUNC_CSS <- norm.func(transposed_df_krakenILUNC, 'CSS')
ILUNC_DEQ <- norm.func(transposed_df_krakenILUNC, 'DeSEQ')
ILUNC_RLE <- norm.func(transposed_df_krakenILUNC, 'RLE+')
ILUNC_RLEpos <- norm.func(transposed_df_krakenILUNC, 'RLE_poscounts')
ILUNC_TMM <- norm.func(transposed_df_krakenILUNC, 'TMM')
ILUNC_logcpm <- norm.func(transposed_df_krakenILUNC, 'logcpm')
ILUNC_rarefy <- norm.func(transposed_df_krakenILUNC, 'rarefy')
ILUNC_CLR <- norm.func(transposed_df_krakenILUNC, 'CLR+')
ILUNC_CLRpos <- norm.func(transposed_df_krakenILUNC, 'CLR_poscounts')
#ILUNC_combat <- norm.func(transposed_df_krakenILUNC, 'combat')


#Adding Stage Labels to Normalization Techniques 
stage_labels <- data.frame(pathologic_stage_label = kraken_metaCOADRNA_Illumina_UNC_clean$pathologic_stage_label,
                           row.names = row.names(kraken_metaCOADRNA_Illumina_UNC_clean))
stage_labels <- t(stage_labels)
#Add Stage Label Function
add_stage_labels <- function(df, stage_labels) {
  common_samples <- intersect(colnames(df), colnames(stage_labels))
  if (length(common_samples) > 0) {
    df <- df[, common_samples]  # Retain only columns present in both dataframes
    df <- rbind(df, stage_labels[, common_samples])
    rownames(df)[nrow(df)] <- "stage_label"
  } else {
    warning("No common samples found between dataframe and stage labels.")
  }
  return(df)
}
#Call Function for Each Normalization Technique
ILUNC_TSS <- add_stage_labels(ILUNC_TSS, stage_labels)
ILUNC_MED <- add_stage_labels(ILUNC_MED, stage_labels)
ILUNC_GMPR <- add_stage_labels(ILUNC_GMPR, stage_labels)
ILUNC_UQ <- add_stage_labels(ILUNC_UQ, stage_labels)
ILUNC_CSS <- add_stage_labels(ILUNC_CSS, stage_labels)
ILUNC_DEQ <- add_stage_labels(ILUNC_DEQ, stage_labels)
ILUNC_RLE <- add_stage_labels(ILUNC_RLE, stage_labels)
ILUNC_RLEpos <- add_stage_labels(ILUNC_RLEpos, stage_labels)
ILUNC_TMM <- add_stage_labels(ILUNC_TMM, stage_labels)
ILUNC_logcpm <- add_stage_labels(ILUNC_logcpm, stage_labels)
ILUNC_rarefy <- add_stage_labels(ILUNC_rarefy, stage_labels)
ILUNC_CLR <- add_stage_labels(ILUNC_CLR, stage_labels)
ILUNC_CLRpos <- add_stage_labels(ILUNC_CLRpos, stage_labels)
#ILUNC_combat <- add_stage_labels(ILUNC_combat, stage_labels)

install.packages("uwot")
library(uwot)
library(ggplot2)
install.packages("forcats")
library(forcats)

combined_data <- cbind(ILUNC_TSS, ILUNC_MED, ILUNC_GMPR, ILUNC_UQ, ILUNC_CSS, ILUNC_DEQ, ILUNC_RLE,
                       ILUNC_RLEpos, ILUNC_TMM, ILUNC_logcpm, ILUNC_rarefy, ILUNC_CLR, ILUNC_CLRpos)
combined_data$stage_label_numeric <- as.numeric(factor(combined_data$stage_label))
combined_data_numeric <- as.matrix(combined_data[-ncol(combined_data)])
combined_data <- t(combined_data_numeric)


umap_result_all <- umap(combined_data[-nrow(combined_data), ], n_neighbors = 15, n_components = 2, metric = "euclidean")

# Plot UMAP with points colored by dataset
ggplot(as.data.frame(umap_result_all$layout), aes(x = V1, y = V2)) +
  geom_point(aes(color = rep(1:ncol(combined_data), each = nrow(combined_data))), size = 2) +
  scale_color_manual(values = rainbow(ncol(combined_data))) +
  ggtitle("UMAP Visualization of All Datasets") +
  theme_minimal()

# Add legend
legend("topright", legend = colnames(combined_data)[-ncol(combined_data)], 
       col = 1:ncol(combined_data) - 1, pch = 20)
