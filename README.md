# 5loclung (5 location lung study)
## A spatially resolved atlas of the human lung characterizes a gland-associated immune niche

<p align="center">

<img src="https://user-images.githubusercontent.com/77395759/226656167-c3c24c2e-0740-47ae-aec9-ab6f0f3ccfc6.png" width=50% height=50%>

</p>

This repository contains custom codes used in the analysis of single cell, nuclei and spatial transcriptomics data from the healthy human lung, now published in [Nature Genetics](https://www.nature.com/articles/s41588-022-01243-4).

Visit our CellxGene browser!: https://www.lungcellatlas.org/

### Code availability

Most of the codes used in manuscript are publicly available packages with specifications written in the methods of the study.

- Code for fGWAS plots and for cell type proportion analysis is available here: https://github.com/natsuhiko/PHM

- Code for marker gene dot plots with mean group expressions and expression of TCR regions were previously published [(Park, J et al. Science 2020)](https://www.science.org/doi/10.1126/science.aay3224) and the code available [here](https://zenodo.org/record/3711134#.ZBnPc-zP2qA) (10.5281/zenodo.3711134)

- Code and data from cell2location analysis of Visium data is available [here](https://github.com/vitkl/adult_lung_mapping/)

- Code for shared TCR clonotype analysis across donors and locations is in the tile TCR-clonotypes.ipynb.

- Code for explained variability code is explained below and in folders Data, Explained Variability and Plots.

### Data Availability

The processed scRNA-seq, snRNA-seq and Visium ST data are available for browsing and download via our website www.lungcellatlas.org. The dataset (raw data and metadata) is available on the [Human Cell Atlas Data Portal](https://data.humancellatlas.org/explore/projects/957261f7-2bd6-4358-a6ed-24ee080d5cfc) and on the European Nucleotide Archive (ENA) under accession number [PRJEB52292](https://www.ebi.ac.uk/ena/browser/view/PRJEB52292) and BioStudies accession [S-SUBS17](https://www.ebi.ac.uk/biostudies/dsp/studies/S-SUBS17). The Visium data are publicly available on ArrayExpress with the accession number [E-MTAB-11640](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11640). Imaging data can be downloaded from European Bioinformatics Institute (EBI) BioImage Archive under accession number [S-BIAD570](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD570). Additional data were accessed to support analysis and conclusions, which can be accessed through National Centre for Biotechnology Information Gene Expression Omnibus [GSE136831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831), and [GSE134174](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134174) and the HLCA integration, which can be accessed at https://github.com/LungCellAtlas/HLCA.

### Instructions on the analysis for calculating explained varibility by a metadata factor

### How to Execute

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
