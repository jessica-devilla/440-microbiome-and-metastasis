#devtools::install_github("wdl2459/ConQuR")install.packages("doParallel")
#devtools::install_github("wdl2459/ConQuR")
#BiocManager::install("sva")
#BiocManager::install("limma")
#BiocManager::install("metagenomeSeq")
#BiocManager::install("edgeR")
#test push

suppressPackageStartupMessages({
  library(doParallel) 
  library(phyloseq)
  library(DESeq2)
  library(limma)
  library(edgeR)
  library(ConQuR)
  library(sva)
  library(metagenomeSeq)
  library(compositions)
  library(GUniFrac)
})


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
