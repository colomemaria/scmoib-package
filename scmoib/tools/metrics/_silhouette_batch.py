from scib.metrics.silhouette import silhouette_batch as scib_sil_batch
from anndata import AnnData
from typing import Tuple


def silhouette_batch(
        adata: AnnData,
        batch_key: str,
        cell_label: str,
        embed: str = 'X_pca',
        metric: str = 'euclidean',
        scale: bool = True
):
    """
    Calculate Average width slhouette score
    
    Parameters
    ----------
    adata
        Annotated data matrix
    batch_key
        Obs variable containing batch annotation
    cell_label
        Obs variable containing cell type annotation
    embed
        Obsm variable name
    metric
        Distance metric
    scale
        If True, score is scaled between 0 and 1
    Returns
    -------
    Average width slhouette score
    """
    return scib_sil_batch(adata, batch_key, cell_label, embed, metric, scale=scale, return_all=False, verbose=False)
