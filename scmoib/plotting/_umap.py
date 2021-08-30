import scanpy as sc


def umap(adata, color=None):
    sc.pl.umap(adata, color=color, edges=True)