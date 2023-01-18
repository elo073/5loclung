# 5loclung (5 location lung study)


Custom codes accompanying single cell spatial transcriptomic study on the healthy human lung.

This repository contains custom codes used in the analysis of single cell, nuclei and spatial transcriptomics data from the healthy human lung.
The pre-print is available here: https://www.biorxiv.org/content/10.1101/2021.11.26.470108v1
The peer-reviewed publication is available here: https://www.nature.com/articles/s41588-022-01243-4


Most of the codes used in manuscript are publicly available packages with specifications written in the methods of the study.

The code used for fGWAS plots and for cell type proportion analysis is available here: https://github.com/natsuhiko/PHM

The code used for marker gene dot plots with mean group expressions and expression of TCR regions were previously published (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7611066/#, doi: 10.1126/science.aay3224) and the code available here: 10.5281/zenodo.3711134



The code for shared TCR clonotype analysis across donors and locations is in the tile TCR-clonotypes.ipynb.





The code for explained variability code is explained below and in folders Data, Explained Variability and Plots.




## Instructions on the analysis for calculating explained varibility by a metadata factor

## How to Execute

First, clone the repository
```bash
$ git clone https://github.com/elo073/5loclung.git
```

Next, access the  data portal (<https://5locationslung.cellgeni.sanger.ac.uk/cellxgene.html>) and download the H5AD object under "All data". Save it in 5loclung/Data

Finaly, run the following commands:

```bash
# Access the script's folder
$ cd 5loclung/Explained\ Variability/ 

# Write count matrices
$ python convert_h5ad.py
$ python convert_h5ad_smg.py

# Execute scripts for explained variability
$ Rscript run.R
$ Rscript run_smg_sc.R
$ Rscript tun_smg_sn.R

```

The plots will be saved in the 'Plots' folder
