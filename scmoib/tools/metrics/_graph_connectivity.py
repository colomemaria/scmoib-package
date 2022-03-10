# This original version of this code was written for the scIB project
# For more information see: https://github.com/theislab/scib
# Paper to cite for this code : https://www.nature.com/articles/s41592-021-01336-8
# M. D. Luecken, M. Bu ̈ttner, K. Chaichoompu, A. Danese, M. Interlandi, M. F. Mueller, D. C. Strobl, L. Zappia,
# M. Dugas, M. Colome ́-Tatche ́, and F. J. Theis. Benchmarking atlas-level data integration in single-cell genomics.
# https://www.nature.com/articles/s41592-021-01336-8

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components
from anndata import AnnData


def graph_connectivity(adata: AnnData, label_key: str) -> float:
    """"
    Metric that quantifies how connected the subgraph corresponding to each batch cluster is.

    Parameters
    ----------
    adata
        Annotated data matrix.
    label_key
        Batch cluster key.
    Returns
    -------
    Mean graph connectivity score.
    """
    if 'connectivities' not in adata.obsp:
        raise KeyError('Please compute the neighborhood graph before running this '
                       'function!')

    adata.obs[label_key] = adata.obs[label_key].astype('category')
    clust_res = []

    for ct in adata.obs[label_key].cat.categories:
        adata_sub = adata[adata.obs[label_key].isin([ct]), ]
        _, labs = connected_components(adata_sub.obsp['connectivities'], connection='strong')
        tab = pd.value_counts(labs)
        clust_res.append(tab.max() / sum(tab))

    return np.mean(clust_res)
