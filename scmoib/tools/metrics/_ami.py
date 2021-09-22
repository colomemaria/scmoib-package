import episcanpy as epi
from anndata import AnnData


def ami(adata: AnnData, group1: str, group2: str) -> float:
    """
    Computes the AMI score.
    Parameters
    ----------
    adata
        Annotated data matrix.
    group1
        Ground-truth cluster assignments (e.g. cell type labels).
    group2
        "predicted" cluster assignments.
    The function is symmetric, so group1 and group2 can be switched
    Returns
    -------
    The AMI score
    """
    return epi.tl.AMI(adata, group1, group2)
