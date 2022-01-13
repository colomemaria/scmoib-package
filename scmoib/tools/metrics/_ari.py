# This original version of this code was written for the scIB project
# For more information see: https://github.com/theislab/scib
# Paper to cite for this code : https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2
# M. D. Luecken, M. Bu ̈ttner, K. Chaichoompu, A. Danese, M. Interlandi, M. F. Mueller, D. C. Strobl, L. Zappia,
# M. Dugas, M. Colome ́-Tatche ́, and F. J. Theis. Benchmarking atlas-level data integration in single-cell genomics.
# bioRxiv, page 2020.05.22.111161, May 2020. doi: 10.1101/2020.05.22.111161.

import scib
from anndata import AnnData


def ari(adata: AnnData, group1: object, group2: object) -> float:
    """
    Computes the ARI score.
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
    The ARI score
    """

    return scib.metrics.ari(adata, group1, group2)
