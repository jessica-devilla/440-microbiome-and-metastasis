library(GUniFrac)
library(dplyr)
library(ggplot2)

# modify the zicoseq plot function to fix color scale


#' Plot ZicoSeq results
#' The function plots the association strength and direction based on the output from \code{ZicoSeq}.
#' @param ZicoSeq.obj return from function \code{ZicoSeq}.
#' @param pvalue.type character; It could be 'p.adj.fdr','p.raw' or 'p.adj.fwer'.
#' @param cutoff a real value between 0 and 1; cutoff for pvalue.type.
#' @param text.size text size for the plots.
#' @param out.dir character; the directory to save the figures, e.g., \code{getwd()}. Default is NULL. If NULL, figures will not be saved.
#' @param file.name character; the name of the file.
#' @param width the width of the graphics region in inches. See R function \code{ggsave}.
#' @param height the height of the graphics region in inches. See R function \code{ggsave}.
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @import dplyr
#' @import tibble
#' @importFrom ggpubr ggarrange
#' @rdname ZicoSeq.plot
#' @export
#' 


ZicoSeq.plot <- function(ZicoSeq.obj, pvalue.type = c('p.adj.fdr','p.raw','p.adj.fwer'), 
                         cutoff = 0.05, text.size = 10, out.dir = NULL, file.name = 'ZicoSeq.plot.pdf',  width = 10, height = 6,
                         color_scale_limits = c(-50, 20), num_labels=10){
  
  #  if(sum(ZicoSeq.obj[[pvalue.type]]<= cutoff) == 0){
  #    cat(paste0('No taxa within ', pvalue.type, ' <= ',cutoff, '!'))
  #  }
  
  grp.name <-  ZicoSeq.obj$grp.name
  meta.dat <- ZicoSeq.obj$meta.dat
  
  # cat('Significant(',pvalue.type,'<=',cutoff,') association strength between taxa and ',grp.name, ' is visualized!\n' )
  
  if(ZicoSeq.obj$call$feature.dat.type == 'proportion'){
    abundance <- ZicoSeq.obj$feature.dat
    prevalence <- apply(ZicoSeq.obj$feature.dat, 1, function(x) mean(x > 0))
  }
  
  if(ZicoSeq.obj$call$feature.dat.type == 'count'){
    abundance <- t(t(ZicoSeq.obj$feature.dat)/colSums(ZicoSeq.obj$feature.dat))
    prevalence <- apply(ZicoSeq.obj$feature.dat, 1, function(x) mean(x > 0))
  }
  
  if(ZicoSeq.obj$call$feature.dat.type == 'other'){
    abundance <- ZicoSeq.obj$feature.dat
    prevalence <- apply(abundance, 1, function(x) sd(x))
  }
  
  # relative abundance
  abundance <- apply(abundance, 1, function(x) mean(x))
  
  if(!(length(unique(meta.dat[,grp.name])) > 2 & !is.numeric(meta.dat[,grp.name]))) {
    coefs <- sapply(ZicoSeq.obj$coef.list, function(x) x[grep(grp.name,rownames(x)),])
    colnames(coefs) <- paste0('coef_Func',1:length(ZicoSeq.obj$coef.list))
  }  
  
  if (!is.numeric(meta.dat[,grp.name]) & length(unique(meta.dat[,grp.name])) == 2) {
    level2 <- rownames(ZicoSeq.obj$coef.list[[1]])[grep(paste0('^',grp.name),rownames(ZicoSeq.obj$coef.list[[1]]))]
    level2 <- gsub(grp.name, '', level2)
    base <- setdiff(unique(meta.dat[,grp.name]), level2)
  }
  
  # Break the ties
  ZicoSeq.obj$R2 <- ZicoSeq.obj$R2 + runif(length(ZicoSeq.obj$R2), 0, 1e-10)
  
  R2 <- rowMaxs(ZicoSeq.obj$R2)
  
  ## pvalue = 0??? -log10(p)
  if(length(unique(meta.dat[,grp.name])) > 2 & !is.numeric(meta.dat[,grp.name])){
    plot.data <- data.frame(pvals = ZicoSeq.obj[[pvalue.type]], 
                            prevalence = prevalence, abundance = abundance, 
                            R2 = R2, 
                            taxa = rownames(ZicoSeq.obj$R2))
    plot.data[plot.data$pvals > cutoff, 'taxa'] <- ''
    
    # Sort data based on p-values
    plot.data <- plot.data[order(plot.data$pvals), ]
    
    # Subset to include only the data with the 10 lowest p-values
    plot.data_top <- plot.data[1:num_labels, ]
    
  }else{
    
    signs <- t(sign(coefs))[t(ZicoSeq.obj$R2 == R2)]
    signs[signs == 0] <- 1
    
    plot.data <- data.frame(pvals = ZicoSeq.obj[[pvalue.type]], 
                            prevalence = prevalence, abundance = abundance, 
                            R2 = R2 * signs, 
                            taxa = rownames(ZicoSeq.obj$R2))
    
    plot.data[plot.data$pvals > cutoff, 'taxa'] <- ''
    
    # Sort data based on p-values
    plot.data <- plot.data[order(plot.data$pvals), ]
    
    # Subset to include only the data with the 10 lowest p-values
    plot.data_top <- plot.data[1:num_labels, ]
  }
  
  if(is.numeric(meta.dat[,grp.name])) {
    title.name = paste0('Association between ', grp.name, ' and  taxa abundance')
  }
  
  if(!is.numeric(meta.dat[,grp.name]) & length(unique(meta.dat[,grp.name])) > 2){
    title.name = paste0('Association between ', grp.name, ' and  taxa abundance')
  }
  
  if(!is.numeric(meta.dat[,grp.name]) & length(unique(meta.dat[,grp.name])) == 2){
    title.name = paste0('Differential abundance between ', paste0(base, ' (reference) and ', level2))
  }
  
  pvals <- R2 <- taxa <- NULL
  plot.obj <-
    ggplot(plot.data, aes(x = R2, y = -log10(pvals))) +
    geom_point(aes(size = abundance, color = log(prevalence))) + ### CHANGE TO LOG PREVALENCE
    geom_vline(aes(xintercept = 0), color = 'gray', linetype = 'dashed') +
    geom_hline(aes(yintercept = -log10(cutoff)), color = 'gray', linetype = 'dashed') +
    scale_colour_gradient(low = "white", high = "#006D2C") +
    scale_y_continuous(limits = c(0, max(-log10(plot.data$pvals)) * 1.3)) +
    ggrepel::geom_text_repel(data =plot.data_top, aes(label = taxa), max.overlaps = Inf, color = 'black',size=3) +
    labs(x = bquote(R^2), y = paste0('-log10(',pvalue.type,')'),
         color = ifelse(ZicoSeq.obj$call$feature.dat.type == 'other','Standard deviation','Prevalence'),
         size = ifelse(ZicoSeq.obj$call$feature.dat.type == 'other','Mean value','Mean abundance')) +
    theme_bw() +
    theme(axis.text = element_text(color = 'black', size = text.size),
          axis.title = element_text(color = 'black', size = text.size),
          legend.text = element_text(color = 'black', size = text.size),
          legend.title = element_text(color = 'black', size = text.size)) +
    ggtitle(title.name)
  
  
  
  #	plots <- ggarrange(plotlist = plot.list, common.legend = T)
  
  if(!is.null(out.dir)) {
    print(plot.obj)
    ggsave(paste0(out.dir, file.name), width = width, height = height)
  }
  return(plot.obj)
}

run_zicoseq <- function(kraken_data, kraken_meta, name){
  
  
  #subset dataframes for comparison
  kraken_meta_sub <- kraken_meta[kraken_meta$pathologic_stage_label %in% c("Stage I", "Stage IV"), ]
  kraken_subset <- kraken_data[rownames(kraken_data) %in% rownames(kraken_meta_sub), ]
  
  
  # examine dataset
  kraken_mat <- as.matrix(kraken_subset)
  kraken_mat <- as.numeric(kraken_mat)
  print(typeof(kraken_mat))
  zero_count <- sum(kraken_subset == 0) 
  print(kraken_subset[kraken_subset==0])
  nan_count <- sum(is.nan(kraken_mat))
  print(kraken_subset[is.nan(kraken_mat)])
  
  
  #remove columns containing all zeroes from dataframe
  kraken_subset <- kraken_subset[, colSums(kraken_subset != 0) > 0]
  
  #remove columns containing less than 3 true values
  non_zero_counts <- colSums(kraken_subset != 0)
  kraken_subset <- kraken_subset[, non_zero_counts >= 3]
  
  #get only genus labels
  colnames(kraken_subset) <- sub(".*__(.*)$", "\\1", colnames(kraken_subset))
  
  zicoseq_data <- t(kraken_subset)
  #row_names <- 1:nrow(zicoseq_data)
  #rownames(zicoseq_data) <- row_names
  
  zicoObj <- ZicoSeq(meta.dat = kraken_meta_sub,  feature.dat = zicoseq_data, grp.name = 'pathologic_stage_label', 
                     adj.name = NULL, feature.dat.type = "other", prev.filter = 0, mean.abund.filter = 0,  max.abund.filter = 0, min.prop = 0, 
                     # Winsorization to replace outliers
                     is.winsor = TRUE, outlier.pct = 0.001, winsor.end = 'top',
                     # Posterior sampling 
                     is.post.sample = TRUE, post.sample.no = 25, 
                     # Use the identity function
                     link.func = list(function (x) x), stats.combine.func = max,
                     # Permutation-based multiple testing correction
                     perm.no = 99,  strata = NULL, 
                     # Reference-based multiple stage normalization
                     ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                     # Family-wise error rate control
                     is.fwer = FALSE,
                     verbose = TRUE, return.feature.dat = TRUE)
  
  zico_plot <- ZicoSeq.plot(zicoObj, pvalue.type = 'p.adj.fdr', cutoff = 0.05, text.size = 10,
                                out.dir = NULL, width = 15, height = 11)
  
  filename <- paste0("figures/zico_seq/zicoseq_stageivsiv_", name, ".png")
  ggsave(filename, plot = zico_plot)
  
  return(list(zicoObj = zicoObj, zicoseq_data = zicoseq_data))
  
}
