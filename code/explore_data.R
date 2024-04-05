#!/usr/bin/env/Rscript --vanilla


# Clean environment -------------------------------------------------------
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# Print a starting message
cat("Running exploratory analysis...\n")

# load the libraries
suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(dplyr)
  library(corrr)
  library(ggcorrplot)
  library(FactoMineR)
  library(devtools)
  library(ggbiplot)
  library(factoextra)
  library(ggrepel)
})

#### EXPLORATORY ANALYSIS OF THE DATASET
kraken_COAD = readRDS('data/kraken_COAD.RDS')
kraken_metaCOAD  =readRDS('data/kraken_metaCOAD.RDS')


# plot stage labels
stage_hist1 <- ggplot(kraken_metaCOAD, aes(x = pathologic_stage_label)) +
  # geom_bar() +
  labs(x = "Pathologic Stage Label", y = "Count") +
  ggtitle("Histogram of Pathologic Stage Labels")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(stage_hist1)

# create new df which combines stage labels and plots histogram
kraken_meta_filt <-kraken_metaCOAD %>%
  mutate(stage_category = case_when(
    grepl("Stage IV[A-C]*", pathologic_stage_label, ignore.case = TRUE) ~ "Stage IV",
    grepl("Stage III[A-C]*", pathologic_stage_label, ignore.case = TRUE) ~ "Stage III",
    grepl("Stage II[A-C]*", pathologic_stage_label, ignore.case = TRUE) ~ "Stage II",
    grepl("Stage I[A-C]*", pathologic_stage_label, ignore.case = TRUE) ~ "Stage I",
    TRUE ~ "Other"
    
  )
  )

# plot stage labels combined
stage_hist2 <- ggplot(kraken_meta_filt, aes(x = stage_category)) +
  geom_bar() +
  labs(x = "Pathologic Stage Label", y = "Count") +
  ggtitle("Histogram of Pathologic Stage Labels")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(stage_hist2)


# pc ----------------------------------------------------------------------


colnames(kraken_meta_filt)[1] <- "id"
colnames(kraken_COAD)[1] <- "id"

kraken_merge <- merge(kraken_COAD,kraken_meta_filt[,c('id','stage_category')], by='id',all.x=TRUE)
contaminant_columns <- grepl("contaminant", names(kraken_merge), ignore.case = TRUE)

# Subset the dataframe to exclude contaminant columns
kraken_merge <- kraken_merge[, !contaminant_columns]

#merge by species
colnames(kraken_merge) <- sub(".*__(.*)$", "\\1", colnames(kraken_merge))
kraken_pca <- kraken_merge[ , !(names(kraken_merge) %in% c("id", "stage_category"))]
kraken_pca <- scale(kraken_pca)

# METHOD 1
#corr_matrix <- cor(kraken_pca)
#data.pca <- princomp(corr_matrix)
#data.pca$loadings[1:10, 1:2]

# METHOD 2
data.pca <- prcomp(kraken_pca, center = FALSE, scale = FALSE)
data.pca$rotation[1:10, 1:2]

summary(data.pca)

fviz_eig(data.pca, addlabels = TRUE)

## GRAPHS OF INDIVIDUALS
# PCA colored by cos2
fviz_pca_ind(data.pca, geom="point",col.ind="cos2")+
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=0.6)

# Color individuals by group
fviz_pca_ind(data.pca, label="none", 
             col.ind = kraken_merge$stage_category)


##GRAPHS OF VARIABLES
#fviz_pca_var(data.pca, select.var = list(contrib = 20))


g <- fviz_pca_biplot(data.pca,
                     label='var', 
                     repel = TRUE, 
                     select.var = list(contrib = 5),
                     alpha.ind = 0.5,
                     col.ind = "gray",
                     col.var = 'blue',
                     labelsize = 3,
                     arrowsize = 1,
                     alpha.var=0.15,
                     #col.ind = kraken_merge$stage_category
)

print(g)

x <- readRDS("data/kraken_df.R")

# subset dataframe function -----------------------------------------------

subset_and_summarize_by_prefix <- function(df, prefix) {
  
  # Subset by prefix
  colnames(df) <- sub(paste0("^.*", prefix, "__([^\\.]+)\\..*"), "\\1", colnames(df))
  
  # Remove contaminant columns
  contaminant_columns <- grepl("contaminant", names(df), ignore.case = TRUE)
  df <- df[, !contaminant_columns]
  
  # Remove virus columns
  virus_columns <- grepl("^k__", names(df), ignore.case = TRUE)
  df <- df[, !virus_columns]
  
  # Get unique phylum names
  phylas <- unique(sapply(strsplit(names(df), "\\."), `[`, 1))
  
  # Prepare a new dataframe to store sums
  sums_df <- data.frame(id = df[, 1])
  
  # Sum like-named columns across rows
  for(phylum in phylas) {
    pattern <- paste0("^", phylum, "\\.")
    cols <- grep(pattern, names(df), value = TRUE)
    
    # Using apply() to sum across rows for matched columns
    if(length(cols) > 0) {
      row_sums <- rowSums(df[, cols, drop = FALSE], na.rm = TRUE)
      sums_df[[phylum]] <- row_sums  # Storing the sums in the new dataframe
    }
  }
  
  return(sums_df)
}

# pca by phylum -----------------------------------------------------------

kraken_phyla = kraken_COAD

# Assuming kraken_phyla is your dataframe
kraken_phyla_sums <- subset_and_summarize_by_prefix(kraken_phyla, "f")
kraken_phyla_sums <- merge(kraken_phyla_sums,kraken_meta_filt[,c('id','stage_category')], by='id',all.x=TRUE)

kraken_pca <- kraken_phyla_sums[ , !(names(kraken_phyla_sums) %in% c("id", "stage_category"))]

kraken_pca <- scale(kraken_pca)

# Run PCA
data.pca <- prcomp(kraken_pca, center = FALSE, scale = FALSE)
data.pca$rotation[1:10, 1:2]

summary(data.pca)

fviz_eig(data.pca, addlabels = TRUE)

##GRAPHS OF VARIABLES
fviz_pca_var(data.pca, select.var = list(contrib = 20))

# Color individuals by group
fviz_pca_ind(data.pca, label="none", 
             col.ind = kraken_merge$stage_category)


g <- fviz_pca_biplot(data.pca,
                     label='var', 
                     repel = TRUE, 
                     select.var = list(contrib = 5),
                     alpha.ind = 0.5,
                     col.ind = "gray",
                     col.var = 'blue',
                     labelsize = 3,
                     arrowsize = 1,
                     alpha.var=0.15
)
print(g)


