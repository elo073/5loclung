import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as de

adata = sc.read_h5ad("../Data/lung_5loc_sc_sn_cellxgene_030222.h5ad")
adata = adata[adata.obs.Celltypes.isin(["SMG_Duct", "SMG_Mucous", "SMG_Serous"])].copy()
adata = adata[adata.obs.Loc_true.isin(['a_Trachea', 'b_Bronchi.2.3', 'c_Bronchi.4'])].copy()

cells = adata[adata.obs.scsn=="cells"].copy()
nuclei = adata[adata.obs.scsn=="nuclei"].copy()
nuclei = nuclei[nuclei.obs.Donor!="A42"].copy()

to_write_sc = pd.DataFrame(de.csr_matrix.todense(cells.X), index=cells.obs.index, columns=cells.var.index).T
to_write_sn = pd.DataFrame(de.csr_matrix.todense(nuclei.X), index=nuclei.obs.index, columns=nuclei.var.index).T


print(to_write_sc.head())
print(to_write_sn.head())


to_write_sc.to_csv("../Data/smg_sc.csv")
to_write_sn.to_csv("../Data/smg_sn.csv")

