if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EDASeq")


library(TCGAbiolinks)
library(dplyr)
library(DT)
library(EDASeq)


## TEST QUERY TO MAKE SURE IT WORKS

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


### NOW DO OUR PROJECT DATA

query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  access = "open",
  sample.type = "Primary Tumor",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
)


datatable(
  getResults(query), 
  filter = 'top',
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)

# downloads everything into data file (481 files, ~2 GB)
#takes a while - saves to parent folder in chunks then moves to subfolder in data TCGA-COAD
GDCdownload(query = query, directory = "data")

# saves result as a summarized experiment (.rda file is 1.1 GB)
se <- GDCprepare(query, 
                 save=TRUE,
                 save.filename = "data/COAD_geneExp_dataframe.rda",
                 summarizedExperiment = TRUE,
                 directory = 'data')

se <- load(file ="data/COAD_geneExp_dataframe.rda", verbose=FALSE)
se <-data
## get gene Expression values
geneExp <- SummarizedExperiment::assay(se)

# get subtype information
infomation.subtype <- TCGAquery_subtype(tumor = "COAD")
# get clinical data
information.clinical <- GDCquery_clinic(project = "TCGA-COAD",type = "clinical") 

## SUBSET GROUPS STAGE I vs STAGE IV
samples.stage.i <-se$barcode[se$ajcc_pathologic_stage=='Stage I']
samples.stage.iv <-se$barcode[se$ajcc_pathologic_stage=='Stage IV']

##remove nans
samples.stage.i <- samples.stage.i[!is.na(samples.stage.i)]
samples.stage.iv <- samples.stage.iv[!is.na(samples.stage.iv)]


dataPrep <- TCGAanalyze_Preprocessing(
  object = se, 
  cor.cut = 0.6
)   
print(dataPrep)

# takes a few minutes
dataNorm <- TCGAanalyze_Normalization(
  tabDF = dataPrep,
  geneInfo = geneInfoHT,
  method = "gcContent"
)
saveRDS(dataNorm, file = "data/gdc_TCGA_norm.RDS")

dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile", 
  qnt.cut =  0.25
)   


dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,samples.stage.i],
  mat2 = dataFilt[,samples.stage.iv],
  Cond1type = "Stage-I",
  Cond2type = "Stage-IV",
  fdr.cut = 0.01 ,
  logFC.cut = 2,
  method = "glmLRT",
  pipeline = "edgeR"
)
saveRDS(dataDEGs, file = "data/gdc_TCGA_stagei_vs_stageiv_DEGs.RDS")

ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes Stage I Vs Stage IV",
  RegulonList = dataDEGs$gene_name
)  

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  GOBPTab = ansEA$ResBP,
  GOCCTab = ansEA$ResCC,
  GOMFTab = ansEA$ResMF,
  PathTab = ansEA$ResPat,
  nRGTab = dataDEGs$gene_name,
  nBar = 10
)

