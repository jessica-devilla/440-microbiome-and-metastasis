#!/usr/bin/env/Rscript --vanilla


# Clean environment -------------------------------------------------------
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# Print a starting message
cat("Running PCA analysis...\n")

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
kraken_COAD <- readRDS('data/kraken_COAD.RDS')
kraken_metaCOAD  <- readRDS('data/kraken_metaCOAD.RDS')



# plot stage labels
stage_hist1 <- ggplot(kraken_metaCOAD, aes(x = pathologic_stage_label)) +
  geom_bar() +
  labs(x = "Pathologic Stage Label", y = "Count") +
  ggtitle("Histogram of Pathologic Stage Labels")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(stage_hist1)
ggsave(path = "figures", filename = "stage_breakdown_histogram.png", bg='white')

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
ggsave(path = "figures", filename = "stage_breakdown_merge_histogram.png", bg='white')


# pca by genus ----------------------------------------------------------------------

colnames(kraken_meta_filt)[1] <- "id"
colnames(kraken_COAD)[1] <- "id"

kraken_merge <- merge(kraken_COAD,kraken_meta_filt[,c('id','stage_category','sample_type')], by='id',all.x=TRUE)
contaminant_columns <- grepl("contaminant", names(kraken_merge), ignore.case = TRUE)

# Subset the dataframe to exclude contaminant columns
kraken_merge <- kraken_merge[, !contaminant_columns]

# exclude viruses
kraken_merge <- remove_viruses(kraken_merge)

#merge by genus
colnames(kraken_merge) <- sub(".*__(.*)$", "\\1", colnames(kraken_merge))
kraken_pca <- kraken_merge[ , !(names(kraken_merge) %in% c("id", "stage_category",'sample_type'))]
kraken_pca <- scale(kraken_pca)


# Run PCA
data.pca <- prcomp(kraken_pca, center = FALSE, scale = FALSE)
data.pca$rotation[1:10, 1:2]

#summary(data.pca)

fviz_eig(data.pca, addlabels = TRUE)

## GRAPHS OF INDIVIDUALS
# PCA colored by cos2
fviz_pca_ind(data.pca, geom="point",col.ind="cos2")+
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=0.6)

# Color individuals by group
fviz_pca_ind(data.pca, label="none",
             col.ind = kraken_merge$stage_category) +
              theme_minimal()+
             scale_shape_manual(values=c(19,19,19,19,19))
#ggsave(path = "figures", filename = "pca_by_genus_withvirus.png", bg='white')


##GRAPHS OF VARIABLES
fviz_pca_var(data.pca, select.var = list(contrib = 10))


g <- fviz_pca_biplot(data.pca,
                     label='var', 
                     repel = TRUE, 
                     select.var = list(contrib = 10),
                     alpha.ind = 0.5,
                     col.ind = "gray",
                     col.var = 'blue',
                     labelsize = 3,
                     arrowsize = 1,
                     alpha.var=0.15,
                     #col.ind = kraken_merge$stage_category
)

print(g)
ggsave(path = "figures", filename = "biplot_by_genus_novirus.png", bg='white')

x <- readRDS("data/kraken_df.R")

# subset dataframe function -----------------------------------------------

remove_viruses <- function(df) {
  contaminant_columns <- grepl("^k__Viruses", names(df), ignore.case = TRUE)
  df <- df[, !contaminant_columns]
  return(df)
}


subset_and_summarize_by_prefix <- function(df, prefix) {
  
  # Subset by prefix
  colnames(df) <- sub(paste0("^.*", prefix, "__([^\\.]+)\\..*"), "\\1", colnames(df))
  
  # Remove contaminant columns
  contaminant_columns <- grepl("contaminant", names(df), ignore.case = TRUE)
  df <- df[, !contaminant_columns]
  
  # Remove columns that dont include prefix
  misc_columns <- grepl("^k__", names(df), ignore.case = TRUE)
  df <- df[, !misc_columns]
  
  # Get unique group names
  prefixes <- unique(sapply(strsplit(names(df), "\\."), `[`, 1))
  print(prefixes)
  
  # Prepare a new dataframe to store sums
  sums_df <- data.frame(id = df[, 1])
  
  # Sum like-named columns across rows
  for(prefix in prefixes) {
    pattern <- paste0("^", prefix, "\\.")
    cols <- grep(pattern, names(df), value = TRUE)
    
    # Using apply() to sum across rows for matched columns
    if(length(cols) > 0) {
      row_sums <- rowSums(df[, cols, drop = FALSE], na.rm = TRUE)
      sums_df[[prefix]] <- row_sums  # Storing the sums in the new dataframe
    }
  }
  
  return(sums_df)
}

# pca by phylum -----------------------------------------------------------

kraken_phyla = kraken_COAD
kraken_phyla =remove_viruses(kraken_phyla)

# Assuming kraken_phyla is your dataframe
kraken_phyla_sums <- subset_and_summarize_by_prefix(kraken_phyla, "p")
kraken_phyla_sums <- merge(kraken_phyla_sums,kraken_meta_filt[,c('id','stage_category')], by='id',all.x=TRUE)

kraken_pca <- kraken_phyla_sums[ , !(names(kraken_phyla_sums) %in% c("id", "stage_category"))]

kraken_pca <- scale(kraken_pca)

# Run PCA
data.pca <- prcomp(kraken_pca, center = FALSE, scale = FALSE)
data.pca$rotation[1:10, 1:2]

#summary(data.pca)

fviz_eig(data.pca, addlabels = TRUE)

##GRAPHS OF VARIABLES
fviz_pca_var(data.pca, select.var = list(contrib = 20))

# Color individuals by group
fviz_pca_ind(data.pca, label="none", 
             col.ind = kraken_merge$stage_category)


g <- fviz_pca_biplot(data.pca,
                     label='var', 
                     repel = TRUE, 
                     select.var = list(contrib = 10),
                     alpha.ind = 0.5,
                     col.ind = "gray",
                     col.var = 'blue',
                     labelsize = 3,
                     arrowsize = 1,
                     alpha.var=0.15
)
print(g)
ggsave(path = "figures", filename = "biplot_by_phylum_novirus.png", bg='white')


# pca by family -----------------------------------------------------------

kraken_fam = kraken_COAD
kraken_fam =remove_viruses(kraken_fam)

# Assuming kraken_phyla is your dataframe
kraken_fam_sums <- subset_and_summarize_by_prefix(kraken_fam, "f")
kraken_fam_sums <- merge(kraken_fam_sums,kraken_meta_filt[,c('id','stage_category')], by='id',all.x=TRUE)

kraken_pca <- kraken_fam_sums[ , !(names(kraken_fam_sums) %in% c("id", "stage_category"))]

kraken_pca <- scale(kraken_pca)

# Run PCA
data.pca <- prcomp(kraken_pca, center = FALSE, scale = FALSE)
data.pca$rotation[1:10, 1:2]

#summary(data.pca)

fviz_eig(data.pca, addlabels = TRUE)

##GRAPHS OF VARIABLES
fviz_pca_var(data.pca, select.var = list(contrib = 20))

# Color individuals by group
fviz_pca_ind(data.pca, label="none", 
             col.ind = kraken_merge$stage_category)


g <- fviz_pca_biplot(data.pca,
                     label='var', 
                     repel = TRUE, 
                     select.var = list(contrib = 10),
                     alpha.ind = 0.5,
                     col.ind = "gray",
                     col.var = 'blue',
                     labelsize = 3,
                     arrowsize = 1,
                     alpha.var=0.15
)
print(g)
ggsave(path = "figures", filename = "biplot_by_family_novirus.png", bg='white')
