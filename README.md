# 440-microbiome-and-metastasis

Jessica Devilla and Grace Goetz

Metastasis is a dynamic and multi-step process by which cancer gains features that allow for invasion of distance sites. Several factors are thought to play a role including tumor immunosuppression mechanisms, metastasis-driving mutations, and phenotypic plasticity, however the mechanistic biology of metastasis is largely unknown. Recently, the intratumoral microbiota has emerged as a driver of cancer progression, and has been utilized in machine learning models to distinguish healthy vs. tumor samples, treatment outcomes, and overall cancer prognosis with impressive predictive power. Ultimately, the use of the cancer-associated microbiome offers a clinical opportunity to predict metastatic spread of cancer and can potentially be modulated to prevent formation of pre-metastatic niches. Given the privileged role of the microbiome in the colon, the high burden of colorectal cancer globally, and association with high rates of metastasis, we chose to examine the impact of the intra-tumoral microbiome on colon adenocarcinoma (COAD). To accomplish this, we used a dataset of genus-level microbial signatures extracted from RNA sequencing and Whole Genome Sequencing data from primary COAD tumors within The Cancer Genome Atlas (TCGA). We found that while there were no global differences in the structure of the microbiome, the abundance of individual strains did vary between stages of COAD. In particular, we found that abundance of Oscillatoria was negatively correlated with progression from stage one to four and was significantly more abundant in stage one tumors. This relationship was conserved across datasets derived from different locations as well as datasets processed using different normalization methods. Taken together with previous studies which suggest Oscillatoria may play a role in  epithelial-to-mesenchymal transition, these results suggest that Oscillatoria may be a strong candidate for diagnostic or therapeutic intervention of metastasis. 

## Data
Poore GD, Kopylova E, Zhu Q, et al. Microbiome analyses of blood and tissues suggest cancer diagnostic approach. Nature. 2020;579(7800):567-574. doi:10.1038/s41586-020-2095-1
- https://ftp.microbio.me/pub/cancer_microbiome_analysis/TCGA/ 
- https://github.com/biocore/tcga

The data contained in the data folder of this repo represents the normalized microbial abudance measurements for a range of cancer types (Kraken-TCGA-Voom-SNM-Plate-Center-Filtering-Data.csv) and corresponding metadata (Metadata-TCGA-Kraken-17625-Samples.csv). This data has been run through a pipeline called Kraken to allow for microbial taxonomy identification and Voom-SNM, a normalization and decontamination process. The result is a csv with values for microbial abundance for each species


## Folder Structure

```
440-microbiome-and-metastasis/
|__ README.md							
|__ code/				
|__ data/						
|__ figures/
```
## Installation

### Prerequisites

This project assumes use of a Unix command line shell and R. Before running, ensure that the Rscript executable is added to environmental PATH variable in order to allow access to Rscript from the command line. For Windows, the default installation location for R is C:\Program Files\R\R-4.3.2\bin. Package dependencies can be found in requirements.txt. Necessary packages will be installed and loaded within in individual scripts.

### Running the code

Clone the github repository into desired directory using Git Bash. 

```bash
git clone https://github.com/jessica-devilla/440-microbiome-and-metastasis.git
```

After cloning the repository, you may need to make the R script executable. Run the following command in the terminal from the root of the cloned repository:
```bash
chmod +x ./code/import_data.R
```

All code to generate figures can be found in the code folder. Each script can be run from bash shell as follows:

```bash
Rscript code/import_data.R
```

import_data.R should be run locally before running additional analysis. This script reads data into R and saves R data files for relevants datasets to the local data folder. They are not included in the repo due to storage requirements.  

The relevant code for generating paper figures are as follows: 

Figure 1: 
```bash
Rscript code/pca_by_phylum_family_genus.R
Rscript code/clusering_bacterial_abundance.R
Rscript code/total_abundances_by_stage.R
Rscript code/phyloseq_alpha_diversity.R
```

Figure 2:
```bash

Rscript code/zicoseq.R
```

Figure 3:
```bash

Rscript code/voom_snm_analysis_figs.R
```

Figure 4:
```bash
Rscript code/phyloseq_by_normalization.R
```

Figure 5:
```bash

```
