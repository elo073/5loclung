# Log CPM matrix (gene x cell)
load("lung.Rbin")
cpm <- data.frame(cpm[,-1], row.names = cpm[,1])

# select genes expressed in at least 5% of all cells
percent <- dim(cpm)[2]*0.95
cpm <- cpm[rowSums(cpm == 0) <= percent, ]

# meta data file
mdata=read.table("Meta_all.csv",header=T, sep=",", row.names=1)

# Logaritmize Total counts and Number of Genes
mdata$Total.counts=log(mdata$Total.counts)
mdata$Number.of.genes=log(mdata$Number.of.genes)

# Linear mixed model
source("data/Sanger/LFLMM.LMMBF.R")

# barplot for the lmm output
source("data/Sanger/Barplot.LFLMM.R")

# fitting lmm
res = LFLMM(cpm, mdata)

# visualising variance explained (%)
Barplot(res, "plot.png", "Cells")

