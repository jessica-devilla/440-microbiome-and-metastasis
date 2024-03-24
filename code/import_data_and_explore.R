
# Load required libraries -------------------------------------------------

## MAKE SURE ALL REQUIREMETS ARE MET AND LOAD LIBRARIES

# Check if the packages is already installed
if (!requireNamespace("readr", quietly = TRUE)) {
  # If not installed, install it
  install.packages("readr")
}
# Check if the package is already installed
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  # If not installed, install it
  install.packages("ggplot2")
}

# Check if the package is already installed
if (!requireNamespace("dplyr", quietly = TRUE)) {
  # If not installed, install it
  install.packages("dplyr")
}

# Check if the package is already installed
if (!requireNamespace("corrr", quietly = TRUE)) {
  # If not installed, install it
  install.packages("corrr")
}

# Check if the package is already installed
if (!requireNamespace("ggcorrplot", quietly = TRUE)) {
  # If not installed, install it
  install.packages("ggcorrplot")
}
# Check if the package is already installed
if (!requireNamespace("FactoMineR", quietly = TRUE)) {
  # If not installed, install it
  install.packages("FactoMineR")
}
# Check if the package is already installed
if (!requireNamespace("ggbiplot", quietly = TRUE)) {
  # If not installed, install it
  install.packages("ggbiplot")
}


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
})


# import data -------------------------------------------------------------

##### IMPORT DATA AND FORMAT

### TEST IMPORT
# read csv files from data folder
#kraken_meta_path <- "C:\Users\jessn\Documents\MIT_G1\20_440\pset6\data\Kraken-TCGA-Voom-SNM-Plate-Center-Filtering-Data.csv"
# set path to csv files via directory
kraken_path <- "data/Kraken-TCGA-Voom-SNM-Plate-Center-Filtering-Data.csv"

#note these import methods do not import the data as in the csv
kraken_data <- read.csv(kraken_path)
#import by directly calling git hub url
kraken_data_git <- read.csv("https://media.githubusercontent.com/media/jessica-devilla/JD_20_440_pset6/main/data/Kraken-TCGA-Voom-SNM-Plate-Center-Filtering-Data.csv")

# import kraken data from poore et al
kraken_url <- "https://media.githubusercontent.com/media/jessica-devilla/JD_20_440_pset6/main/data/Kraken-TCGA-Voom-SNM-Plate-Center-Filtering-Data.csv"
#kraken_data <- read_csv(url(kraken_url),col_types = cols(.default = col_character()))
kraken_data <- read_csv(url(kraken_url),show_col_types = FALSE)
kraken_df <- as.data.frame(kraken_data, stringsAsFactors = FALSE)

# import kraken metadata from poore et al
kraken_meta_url <- "https://media.githubusercontent.com/media/jessica-devilla/JD_20_440_pset6/main/data/Metadata-TCGA-Kraken-17625-Samples.csv"
#kraken_metadata <-read_csv(url(kraken_meta_url),col_types = cols(.default = col_character()))
kraken_metadata <-read_csv(url(kraken_meta_url),show_col_types = FALSE)
kraken_metadata_df <- as.data.frame(kraken_metadata, stringsAsFactors = FALSE)

# Subset the dataframe to get only values from COAD patients
kraken_metaCOAD <- subset(kraken_metadata_df, disease_type == "Colon Adenocarcinoma")
#dim(kraken_metaCOAD)
#get the patient ids from COAD metadata
ids <- kraken_metaCOAD[,1]
kraken_COAD <- kraken_df[kraken_df[,1] %in% ids,]
#dim(kraken_COAD) # check to see if dimensions matched

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

#NOT WORKING YET

kraken_norm = scale(kraken_COAD[,-1])
corr_matrix <- cor(kraken_norm)

data.pca <- princomp(corr_matrix)
summary(data.pca)
data.pca$loadings[1:10, 1:2]
fviz_eig(data.pca, addlabels = TRUE)


#pc <- prcomp(kraken_norm, center=TRUE, scale=TRUE)
#attributes(pc)
#summary(pc)

#g <- ggbiplot(pc,
#              obs.scale = 1,
#              var.scale = 1)
#g <- g + scale_color_discrete(name = '')
#g <- g + theme(legend.direction = 'horizontal',
#            legend.position = 'top')
#print(g)
