from scib.metrics.silhouette import silhouette as scib_sil
from anndata import AnnData
from typing import Tuple


def silhouette(
        adata: AnnData,
        batch_key: str,
        embed: str = 'X_pca',
        metric: str = 'euclidean',
        scale: bool = True
):
    """
    
    Parameters
    ----------
    adata
        Annotated data matrix
    batch_key
        Obs variable containing batch annotation
    embed
        Obsm variable name
    metric
        Distance metric
    scale
        If True, score is scaled between 0 and 1
    Returns
    -------
    Silhouette score
    """
    return scib_sil(adata, batch_key, embed, metric, scale)
