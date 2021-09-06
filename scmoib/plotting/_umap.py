import scanpy as sc


def umap_barcodes(adata, color=None):
    """
    instead of plotting the edge of the knn graph what we want is to draw an 
    edge between the matching barcodes.

    @ Anna : DO IT DAMNED
    """
    sc.pl.umap(adata, color=color, edges=True)