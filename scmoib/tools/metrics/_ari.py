# This original version of this code was written for the scIB project
# For more information see: https://github.com/theislab/scib
# Paper to cite for this code : https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2
# M. D. Luecken, M. Bu ̈ttner, K. Chaichoompu, A. Danese, M. Interlandi, M. F. Mueller, D. C. Strobl, L. Zappia,
# M. Dugas, M. Colome ́-Tatche ́, and F. J. Theis. Benchmarking atlas-level data integration in single-cell genomics.
# bioRxiv, page 2020.05.22.111161, May 2020. doi: 10.1101/2020.05.22.111161.

from .utils.utils import checkAdata, checkBatch
import sklearn
import pandas as pd
from anndata import AnnData


def ari(adata: object, group1: object, group2: object) -> float:
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

    checkAdata(adata)

    if isinstance(group1, str):
        checkBatch(group1, adata.obs)
        group1 = adata.obs[group1].tolist()
    elif isinstance(group1, pd.Series):
        group1 = group1.tolist()

    if isinstance(group2, str):
        checkBatch(group2, adata.obs)
        group2 = adata.obs[group2].tolist()
    elif isinstance(group2, pd.Series):
        group2 = group2.tolist()

    if len(group1) != len(group2):
        raise ValueError(f'different lengths in group1 ({len(group1)}) and group2 ({len(group2)})')

    return sklearn.metrics.cluster.adjusted_rand_score(group1, group2)
