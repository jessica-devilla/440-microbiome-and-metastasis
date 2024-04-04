# 440-microbiome-and-metastasis

Jessica Devilla and Grace Goetz

Data sources
Poore GD, Kopylova E, Zhu Q, et al. Microbiome analyses of blood and tissues suggest cancer diagnostic approach. Nature. 2020;579(7800):567-574. doi:10.1038/s41586-020-2095-1
- https://ftp.microbio.me/pub/cancer_microbiome_analysis/TCGA/ 
- https://github.com/biocore/tcga

Hermida LC, Gertz EM, Ruppin E. Predicting cancer prognosis and drug response from the tumor microbiome. Nat Commun. 2022;13(1):2896. doi:10.1038/s41467-022-30512-3
- https://zenodo.org/records/6471321
- https://github.com/ruppinlab/tcga-microbiome-prediction/tree/v1.2 

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
