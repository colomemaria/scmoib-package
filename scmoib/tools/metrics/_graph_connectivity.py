## This original version of this code was written for the scIB project
# For more information see: https://github.com/theislab/scib
# Paper to cite for this code : https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2
# M. D. Luecken, M. Bu ̈ttner, K. Chaichoompu, A. Danese, M. Interlandi, M. F. Mueller, D. C. Strobl, L. Zappia, M. Dugas, M. Colome ́-Tatche ́, and F. J. Theis. Benchmarking atlas-level data integration in single-cell genomics. bioRxiv, page 2020.05.22.111161, May 2020. doi: 10.1101/2020.05.22.111161.

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components


def graph_connectivity(adata_post, label_key):
    """"
    Metric that quantifies how connected the subgraph corresponding to each batch cluster is.
    """
    if 'connectivities' not in adata_post.obsp:
        raise KeyError('Please compute the neighborhood graph before running this '
                       'function!')

    adata_post.obs[label_key] = adata_post.obs[label_key].astype('category')
    clust_res = []

    for ct in adata_post.obs[label_key].cat.categories:
        adata_post_sub = adata_post[adata_post.obs[label_key].isin([ct]),]
        _, labs = connected_components(adata_post_sub.obsp['connectivities'], connection='strong')
        tab = pd.value_counts(labs)
        clust_res.append(tab[0] / sum(tab))

    return np.mean(clust_res)
