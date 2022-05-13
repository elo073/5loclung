# Linear mixed mode
source("LFLMM.LMMBF.R")

# barplot for the lmm output
source("Barplot.LFLMM.R")


# Log CPM matrix (gene x cell)
load("smg_sc.Rbin")

######select genes expressed in at least 10% of all cells
check_genes = c("CCL28", "CCL20", "CCL2", "TNFSF13", "IL6", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DPB1", "HLA-DQB1", "HLA-DPA1", "HLA-DMA", "HLA-DMB", "HLA-DQA1", "HLA-DPB1", "CD40", "PIGR")

percent <- dim(cpm)[2]*0.95
cpm_filt <- cpm[(rowSums(cpm == 0) <= percent) | (rownames(cpm) %in% check_genes), ]

for (i in check_genes){
	print(paste0(i, " ", i %in% rownames(cpm_filt)))
}

###### Select predetermined genes
cpm_genes <-cpm[rownames(cpm) %in% check_genes, ]


###### meta data file
mdata=read.table("Meta_smg_sc.csv",header=T, sep=",", row.names=1)

# Logaritmize Total counts and Number of Genes
mdata$Total.counts=log(mdata$Total.counts)
mdata$Number.of.genes=log(mdata$Number.of.genes)

####### LOCATION AS LINEAR #########
mdata$Location[mdata$Location == "a_Trachea"] <- 1
mdata$Location[mdata$Location == "b_Bronchi.2.3"] <- 2
mdata$Location[mdata$Location == "c_Bronchi.4"] <- 3

# fitting lmm
res = LFLMM(cpm_filt, mdata)
Barplot(res, "smg_sc_lin.png", "SMG Cells")


###### LOCATION AS BINARY #######
mdata$Location[mdata$Location == 2] <- 0
mdata$Location[mdata$Location == 3] <- 0

# fitting lmm
res = LFLMM(cpm_filt, mdata)
Barplot(res, "smg_sc_bin.png", "SMG Cells")
