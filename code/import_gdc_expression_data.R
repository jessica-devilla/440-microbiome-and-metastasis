if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")


library(TCGAbiolinks)
library(dplyr)
library(DT)


## TEST QUERY

# Gene expression aligned against hg38
query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01")
)
GDCdownload(query = query, directory = "data")
data <- GDCprepare(query = query, directory ='data')

# saves result as a dataframe
df <- GDCprepare(query, 
                 save=TRUE,
                 save.filename = "GBM_geneExp_dataframe.rda",
                 summarizedExperiment = FALSE, 
                 directory ='data')

# saves result as a summarized experiment
se <- GDCprepare(query, 
                 save=TRUE,
                 save.filename = "GBM_geneExp_dataframe.rda",
                 summarizedExperiment = TRUE,
                 directory = 'data')
## get gene Expression values
geneExp <- SummarizedExperiment::assay(se)

