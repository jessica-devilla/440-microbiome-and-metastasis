

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


# import data -------------------------------------------------------------

##### IMPORT DATA AND FORMAT


# import kraken data from poore et al
kraken_url <- "https://media.githubusercontent.com/media/jessica-devilla/JD_20_440_pset6/main/data/Kraken-TCGA-Voom-SNM-Plate-Center-Filtering-Data.csv"
kraken_data <- read_csv(url(kraken_url),show_col_types = FALSE)
kraken_df <- as.data.frame(kraken_data, stringsAsFactors = FALSE)

# import kraken metadata from poore et al
kraken_meta_url <- "https://media.githubusercontent.com/media/jessica-devilla/JD_20_440_pset6/main/data/Metadata-TCGA-Kraken-17625-Samples.csv"
kraken_metadata <-read_csv(url(kraken_meta_url),show_col_types = FALSE)
kraken_metadata_df <- as.data.frame(kraken_metadata, stringsAsFactors = FALSE)

# Save the dataframes as an R file
saveRDS(kraken_df, file = "data/kraken_df.R")
saveRDS(kraken_df, file = "data/kraken_metadatadf.R")

# Subset the dataframe to get only values from COAD patients
kraken_metaCOAD <- subset(kraken_metadata_df, disease_type == "Colon Adenocarcinoma")
#dim(kraken_metaCOAD)
#get the patient ids from COAD metadata
ids <- kraken_metaCOAD[,1]
kraken_COAD <- kraken_df[kraken_df[,1] %in% ids,]
#dim(kraken_COAD) # check to see if dimensions matched

# Save the dataframes as an R file
saveRDS(kraken_COAD, file = "data/kraken_COAD.R")
saveRDS(kraken_metaCOAD, file = "data/kraken_metaCOAD.R")

#### EXPLORATORY ANALYSIS OF THE DATASET

# plot stage labels
stage_hist1 <- ggplot(kraken_metaCOAD, aes(x = pathologic_stage_label)) +
  geom_bar() +
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
