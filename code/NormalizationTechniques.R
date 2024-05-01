#Rarefy Normalization 
install.packages("phyloseq")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Use BiocManager to install ConQuR
BiocManager::install("ConQuR")
BiocManager::install("sva")
BiocManager::install("limma")
BiocManager::install("pamr")


library(phyloseq)
ps <- otu_table(as.matrix(kraken_raw_COADRNA_IlluminaGA_UNC_clean), taxa_are_rows = FALSE)

# Rarefy the data to a specified depth
rarefied_ps <- rarefy_even_depth(ps, sample.size = 10000)
rarefied_COAD_IlluminaUNC <- as.data.frame(otu_table(rarefied_ps))


#DESeq2 Normalization
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

transposed_df_krakenILUNC <- t(kraken_raw_COADRNA_IlluminaGA_UNC_clean)
sample_ids <- colnames(transposed_df_krakenILUNC)
colData <- data.frame(sampleID = sample_ids)

# Create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = transposed_df,
                              colData = colData,
                              design = ~1)
# Estimate size factors
dds <- estimateSizeFactors(dds)
# Retrieve normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
DESeq_COAD_IlluminUNC <- as.data.frame(normalized_counts)






norm.func <- function(p1,norm_method){
  # remove the all zero genes
  p1 <- p1[rowSums(p1)>0,]
  # remove samples with only 1 non-zero features
  p1 <- p1[,colSums(p1!=0)>1]
  # make names for genes
  rownames(p1) <- make.names(rownames(p1))
  
  
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
    rarefied_ps <- rarefy_even_depth(ps, sample.size = 10)
    norm_p1 <- as.data.frame(otu_table(rarefied_ps))
    return(norm_p1)
  }
  
  
  # CLR+, for samples
  # based on TSS normalized data, add a pseudo count 0.65*minimum to zero values
  if(norm_method=="CLR+"){
    require(compositions)
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    # replaced all the 0 abundances with 0.65 times minimum non-zero abundance
    norm_p1[norm_p1==0] <- min(norm_p1[norm_p1!=0])*0.65
    # clr transformation
    trans_p1 <- as.data.frame(apply(norm_p1,2,function(x) clr(x)))
    return(trans_p1)
  }
  

  if(norm_method=="BMC"){
    require(pamr)
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    # Replace all the 0 abundances with 0.65 times the minimum non-zero abundance
    norm_p1[norm_p1==0] <- min(norm_p1[norm_p1!=0])*0.65
    # Log transformation
    trans_p1 <- log(norm_p1)
    # Batch factor
    batch_factor <- factor(rep(1, ncol(trans_p1)))
    # BMC correction
    correct_merged <- as.data.frame(pamr.batchadjust(list(x=as.matrix(trans_p1), batchlabels=batch_factor))$x)
  }

  if(norm_method=="limma"){
    require(limma)
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    # replaced all the 0 abundances with 0.65 times minimum non-zero abundance
    norm_p1[norm_p1==0] <- min(norm_p1[norm_p1!=0])*0.65
    # log transformation
    trans_p1 <- log(norm_p1)
    #limma
    batch_factor <- factor(rep(1, ncol(trans_p1)))
    # Limma normalization
    correct_merged <- as.data.frame(removeBatchEffect(trans_p1, batch=batch_factor))
  }
  
  if (norm_method == "combat") {
    require(sva)
    
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1, 2, function(x) x / sum(x)))
    
    # Replace all the 0 abundances with 0.65 times the minimum non-zero abundance
    norm_p1[norm_p1 == 0] <- min(norm_p1[norm_p1 != 0]) * 0.65
    
    # Log transformation
    trans_p1 <- log(norm_p1)
    
    # Check for genes with no variation within a single batch
    zero_var_genes <- apply(trans_p1, 1, function(row) all(row == 0))
    trans_p1 <- trans_p1[!zero_var_genes, ]
    batch_factor <- factor(rep(1, ncol(trans_p1)))
    correct_combat <- ComBat(trans_p1, batch = batch_factor, ref.batch = 1)
    correct_combat_df <- as.data.frame(correct_combat)
  }
  
  
  if(norm_method=="conqur"){
    require(ConQuR)
    # Batch factor
    batch_factor <- factor(rep(1, ncol(p1)))
    # Covariates (no additional information used)
    covariates <- data.frame(covariate=rep(1, ncol(p1)))
    # ConQuR normalization
    correct_p1 <- ConQuR(tax_tab=as.data.frame(t(p1)), 
                         batchid=batch_factor,
                         batch_ref="1",
                         covariates=covariates,
                         simple_match=TRUE)
    correct_p1 <- as.data.frame(t(correct_p1))
  }
  
}
  
ILUNC_TSS <- norm.func(transposed_df_krakenILUNC, 'TSS')
ILUNC_UQ <- norm.func(transposed_df_krakenILUNC, 'UQ')
ILUNC_CSS <- norm.func(transposed_df_krakenILUNC, 'CSS')
#ILUNC_TMM <- norm.func(transposed_df_krakenILUNC, 'TMM')
ILUNC_logcpm <- norm.func(transposed_df_krakenILUNC, 'logcpm')
ILUNC_rarefy <- norm.func(transposed_df_krakenILUNC, 'rarefy')
ILUNC_CLR <- norm.func(transposed_df_krakenILUNC, 'CLR+')
#ILUNC_BMC <- norm.func(transposed_df_krakenILUNC, 'BMC')
#ILUNC_limma <- norm.func(transposed_df_krakenILUNC, 'limma')
#ILUNC_combat <- norm.func(transposed_df_krakenILUNC, 'combat')
#ILUNC_conqur <- norm.func(transposed_df_krakenILUNC, 'conqur')
