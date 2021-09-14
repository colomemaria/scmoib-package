## This original version of this code was written for the scIB project
# For more information see: https://github.com/theislab/scib
# Paper to cite for this code : https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2
# M. D. Luecken, M. Bu ̈ttner, K. Chaichoompu, A. Danese, M. Interlandi, M. F. Mueller, D. C. Strobl, L. Zappia, M. Dugas, M. Colome ́-Tatche ́, and F. J. Theis. Benchmarking atlas-level data integration in single-cell genomics. bioRxiv, page 2020.05.22.111161, May 2020. doi: 10.1101/2020.05.22.111161.

import sklearn
import pandas as pd
import numpy as np
import scanpy as sc
from ._nmi import nmi


def __opt_louvain(adata, label_key, cluster_key, function=None, resolutions=None,
                  use_rep=None,
                  inplace=True, force=True, verbose=True, **kwargs):
    """
    params:
        label_key: name of column in adata.obs containing biological labels to be
            optimised against
        cluster_key: name of column to be added to adata.obs during clustering. 
            Will be overwritten if exists and `force=True`
        function: function that computes the cost to be optimised over. Must take as
            arguments (adata, group1, group2, **kwargs) and returns a number for maximising
        resolutions: list if resolutions to be optimised over. If `resolutions=None`,
            default resolutions of 20 values ranging between 0.1 and 2 will be used
        use_rep: key of embedding to use only if adata.uns['neighbors'] is not defined,
            otherwise will be ignored
    returns:
        res_max: resolution of maximum score
        score_max: maximum score
        score_all: `pd.DataFrame` containing all scores at resolutions. Can be used to plot the score profile.
        clustering: only if `inplace=False`, return cluster assignment as `pd.Series`
        plot: if `plot=True` plot the score profile over resolution
    """

    if function is None:
        function = nmi

    if cluster_key in adata.obs.columns:
        if force:
            if verbose:
                print(f"Warning: cluster key {cluster_key} already exists " +
                      "in adata.obs and will be overwritten")
        else:
            raise ValueError(f"cluster key {cluster_key} already exists in " +
                             "adata, please remove the key or choose a different name." +
                             "If you want to force overwriting the key, specify `force=True`")

    if resolutions is None:
        n = 20
        resolutions = [2 * x / n for x in range(1, n + 1)]

    score_max = 0
    res_max = resolutions[0]
    clustering = None
    score_all = []

    try:
        adata.uns['neighbors']
    except KeyError:
        if verbose:
            print('computing neigbours for opt_cluster')
        sc.pp.neighbors(adata, use_rep=use_rep)

    for res in resolutions:
        sc.tl.louvain(adata, resolution=res, key_added=cluster_key)
        score = function(adata, label_key, cluster_key, **kwargs)
        score_all.append(score)
        if score_max < score:
            score_max = score
            res_max = res
            clustering = adata.obs[cluster_key]
        del adata.obs[cluster_key]

    if verbose:
        print(f'optimised clustering against {label_key}')
        print(f'optimal cluster resolution: {res_max}')
        print(f'optimal score: {score_max}')

    score_all = pd.DataFrame(zip(resolutions, score_all), columns=('resolution', 'score'))

    if inplace:
        adata.obs[cluster_key] = clustering
        return res_max, score_max, score_all
    else:
        return res_max, score_max, score_all, clustering


def __silhouette(adata, group_key, metric='euclidean', embed='X_pca', scale=True):
    """
    wrapper for sklearn silhouette function values range from [-1, 1] with 1 being an ideal fit, 0 indicating overlapping clusters and -1 indicating misclassified cells
    """
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f'{embed} not in obsm')
    asw = sklearn.metrics.silhouette_score(adata.obsm[embed], adata.obs[group_key], metric=metric)
    if scale:
        asw = (asw + 1) / 2
    return asw


def __silhouette_batch(adata, batch_key, group_key, metric='euclidean',
                       embed='X_pca', verbose=True, scale=True):
    """
    Silhouette score of batch labels subsetted for each group.
    params:
        batch_key: batches to be compared against
        group_key: group labels to be subsetted by e.g. cell type
        metric: see sklearn silhouette score
        embed: name of column in adata.obsm
    returns:
        all scores: absolute silhouette scores per group label
        group means: if `mean=True`
    """
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f'{embed} not in obsm')

    sil_all = pd.DataFrame(columns=['group', 'silhouette_score'])

    for group in adata.obs[group_key].unique():
        adata_group = adata[adata.obs[group_key] == group]
        if adata_group.obs[batch_key].nunique() == 1:
            continue
        sil_per_group = sklearn.metrics.silhouette_samples(adata_group.obsm[embed], adata_group.obs[batch_key],
                                                           metric=metric)
        # take only absolute value
        sil_per_group = [abs(i) for i in sil_per_group]
        if scale:
            # scale s.t. highest number is optimal
            sil_per_group = [1 - i for i in sil_per_group]
        d = pd.DataFrame({'group': [group] * len(sil_per_group), 'silhouette_score': sil_per_group})
        sil_all = sil_all.append(d)
    sil_all = sil_all.reset_index(drop=True)
    sil_means = sil_all.groupby('group').mean()

    if verbose:
        print(f'mean silhouette per cell: {sil_means}')
    return sil_all, sil_means


def __isolated_labels(adata, label_key, batch_key, cluster_key="iso_cluster",
                      cluster=True, n=None, all_=False, verbose=True):
    """
    score how well labels of isolated labels are distiguished in the dataset by
        1. clustering-based approach
        2. silhouette score
    params:
        cluster: if True, use clustering approach, otherwise use silhouette score approach
        n: max number of batches per label for label to be considered as isolated.
            if n is integer, consider labels that are present for n batches as isolated
            if n=None, consider minimum number of batches that labels are present in
        all_: return scores for all isolated labels instead of aggregated mean
    return:
        by default, mean of scores for each isolated label
        retrieve dictionary of scores for each label if `all_` is specified
    """

    scores = {}
    isolated_labels = __get_isolated_labels(adata, label_key, batch_key, cluster_key, n=n, verbose=verbose)
    for label in isolated_labels:
        score = __score_isolated_label(adata, label_key, cluster_key, label, cluster=cluster, verbose=verbose)
        scores[label] = score

    if all_:
        return scores
    return np.mean(list(scores.values()))


def __get_isolated_labels(adata, label_key, batch_key, cluster_key, n, verbose):
    """
    get labels that are considered isolated by the number of batches
    """

    tmp = adata.obs[[label_key, batch_key]].drop_duplicates()
    batch_per_lab = tmp.groupby(label_key).agg({batch_key: "count"})

    # threshold for determining when label is considered isolated
    if n is None:
        n = batch_per_lab.min().tolist()[0]

    if verbose:
        print(f"isolated labels: no more than {n} batches per label")

    labels = batch_per_lab[batch_per_lab[batch_key] <= n].index.tolist()
    if len(labels) == 0 and verbose:
        print(f"no isolated labels with less than {n} batches")
    return labels


def __score_isolated_label(adata, label_key, cluster_key, label, cluster=True, verbose=False, **kwargs):
    """
    compute label score for a single label
    params:
        cluster: if True, use clustering approach, otherwise use silhouette score approach
    """
    adata_tmp = adata.copy()

    def max_label_per_batch(adata, label_key, cluster_key, label, argmax=False):
        """cluster optimizing over cluster with largest number of isolated label per batch"""
        sub = adata.obs[adata.obs[label_key] == label].copy()
        label_counts = sub[cluster_key].value_counts()
        if argmax:
            return label_counts.index[label_counts.argmax()]
        return label_counts.max()

    def max_f1(adata, label_key, cluster_key, label, argmax=False):
        """cluster optimizing over largest F1 score of isolated label"""
        obs = adata.obs
        max_cluster = None
        max_f1 = 0
        for cluster in obs[cluster_key].unique():
            y_pred = obs[cluster_key] == cluster
            y_true = obs[label_key] == label
            f1 = sklearn.metrics.f1_score(y_pred, y_true)
            if f1 > max_f1:
                max_f1 = f1
                max_cluster = cluster
        if argmax:
            return max_cluster
        return max_f1

    if cluster:
        __opt_louvain(adata_tmp, label_key, cluster_key, function=max_f1, label=label, verbose=False, inplace=True)
        score = max_f1(adata_tmp, label_key, cluster_key, label, argmax=False)
    else:
        adata_tmp.obs['group'] = adata_tmp.obs[label_key] == label
        score = __silhouette(adata_tmp, group_key='group', **kwargs)

    del adata_tmp

    if verbose:
        print(f"{label}: {score}")

    return score


def silhouette(adata, batch_key='orig.ident', cell_label='paper.cell.type', embed='X_pca'):
    # global silhouette coefficient
    sil_global = __silhouette(adata, group_key=cell_label, embed=embed, metric='euclidean')
    # silhouette coefficient per batch
    _, sil_clus = __silhouette_batch(adata, batch_key=batch_key, group_key=cell_label,
                                     embed=embed, metric='euclidean', verbose=False)
    sil_clus = sil_clus['silhouette_score'].mean()
    il_score_sil = __isolated_labels(adata, label_key=cell_label, batch_key=batch_key,
                                     cluster=False, n=1, verbose=False)
    il_score_clus = __isolated_labels(adata, label_key=cell_label, batch_key=batch_key,
                                      cluster=True, n=1, verbose=False)

    return sil_global, sil_clus, il_score_clus, il_score_sil
