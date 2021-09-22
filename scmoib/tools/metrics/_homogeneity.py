import episcanpy as epi
from anndata import AnnData


def homogeneity(adata: AnnData, group1: str, group2: str) -> float:
    """
    Computes the homogeneity score.
    Parameters
    ----------
    adata
        Annotated data matrix.
    group1
        #TODO
    group2
        #TODO

    Returns
    -------
    The homogeneity score
    """
    return epi.tl.homogeneity(adata, group1, group2)
