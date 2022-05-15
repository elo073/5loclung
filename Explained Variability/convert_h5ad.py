import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as de

adata = sc.read_h5ad("../Data/lung_5loc_sc_sn_cellxgene_030222.h5ad")
cells = adata[adata.obs.scsn=="cells"].copy()

to_write = pd.DataFrame(de.csr_matrix.todense(cells.X), index=cells.obs.index, columns=cells.var.index).T


print(to_write.head())


to_write.to_csv("../Data/cells.csv")

