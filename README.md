# 5loclung
Custom codes accompanying single cell spatial transcriptomic study on the healthy human lung

This repository contains custom codes used in the analysis of single cell, nuclei and spatial transcriptomics data from the healthy human lung.
The pre-print is available here: https://www.biorxiv.org/content/10.1101/2021.11.26.470108v1

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
