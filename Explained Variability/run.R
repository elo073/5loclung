# Log CPM matrix (gene x cell)
# Execute only the first time
cpm=read.table("../Data/cells.csv", header=T, sep=",")
cpm <- data.frame(cpm[,-1], row.names = cpm[,1])

# select genes expressed in at least 5% of all cells
percent <- dim(cpm)[2]*0.95
cpm <- cpm[rowSums(cpm == 0) <= percent, ]

# meta data file
mdata=read.table("../Data/Meta_all.csv",header=T, sep=",", row.names=1)

# Logaritmize Total counts and Number of Genes
mdata$Total.counts=log(mdata$Total.counts)
mdata$Number.of.genes=log(mdata$Number.of.genes)

# Linear mixed model
source("LFLMM.LMMBF.R")

# barplot for the lmm output
source("Barplot.LFLMM.R")

# fitting lmm
res = LFLMM(cpm, mdata)

# visualising variance explained (%)
Barplot(res, "../Plots/All_sc.png", "Cells")

