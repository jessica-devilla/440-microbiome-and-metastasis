# 440-microbiome-and-metastasis

Jessica Devilla and Grace Goetz

Data sources
Poore GD, Kopylova E, Zhu Q, et al. Microbiome analyses of blood and tissues suggest cancer diagnostic approach. Nature. 2020;579(7800):567-574. doi:10.1038/s41586-020-2095-1
- https://ftp.microbio.me/pub/cancer_microbiome_analysis/TCGA/ 
- https://github.com/biocore/tcga

Hermida LC, Gertz EM, Ruppin E. Predicting cancer prognosis and drug response from the tumor microbiome. Nat Commun. 2022;13(1):2896. doi:10.1038/s41467-022-30512-3
- https://zenodo.org/records/6471321
- https://github.com/ruppinlab/tcga-microbiome-prediction/tree/v1.2

Citation for TCGA Biolinks
Colaprico, Antonio, et al. “TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data.” Nucleic acids research 44.8 (2015): e71-e71.
Silva, Tiago C., et al. “TCGA Workflow: Analyze cancer genomics and epigenomics data using Bioconductor packages.” F1000Research 5 (2016). (https://f1000research.com/articles/5-1542/v2)
Mounir, Mohamed, et al. “New functionalities in the TCGAbiolinks package for the study and integration of cancer data from GDC and GTEx.” PLoS computational biology 15.3 (2019): e1006701. (https://doi.org/10.1371/journal.pcbi.1006701)

Clone the github repository into desired directory using Git Bash. 

```bash
git clone https://github.com/jessica-devilla/440-microbiome-and-metastasis.git
```

After cloning the repository, you may need to make the R script executable. Run the following command in the terminal from the root of the cloned repository:
```bash
chmod +x ./code/impot_data.R
```

To import Poore et al data and save R data files to your local machine, run:

```bash
Rscript code/import_data.R
```
