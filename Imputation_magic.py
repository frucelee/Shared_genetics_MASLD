##imputed gene expression using magic
import scanpy as sc
adata=sc.read_h5ad('Hepatocyte_control_MASLD.h5ad')
import scanpy as sc
import scanpy.external as sce
top_genes_indices=["SUOX","CTCF"]
adata_magic = sce.pp.magic(adata, name_list=["SUOX","CTCF"], knn=5)
adata_magic.shape
import numpy as np
import pandas as pd
import anndata as ad
mat = adata_magic.X
mat = mat.toarray() if hasattr(mat, "toarray") else np.asarray(mat)
df = pd.DataFrame(
mat,
index=adata_magic.obs_names,
columns=adata_magic.var_names
)
df_out = df.copy()
df_out["group"] = adata_magic.obs["group"].values
df_out.to_csv("magic_with_meta.tsv", sep="\t")
