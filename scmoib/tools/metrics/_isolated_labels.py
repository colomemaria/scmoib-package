import pandas as pd
from sklearn.metrics import f1_score
from .silhouette import silhouette
from anndata import AnnData


def isolated_labels(
        adata,
        label_key,
        batch_key,
        embed,
        cluster=True,
        iso_threshold=None,
        return_all=False,
        verbose=True
):
    """
    Score how well labels of isolated labels are distiguished in the dataset by either
        1. clustering-based approach F1 score, or
        2. average-width silhouette score (ASW) on isolated label vs all other labels
    :param adata: anndata object
    :param label_key: column in adata.obs
    :param batch_key: column in adata.obs
    :param cluster: if True, use clustering approach, otherwise use silhouette score approach
    :param embed: key in adata.obsm used for silhouette score if cluster=False, or
        as representation for clustering (if neighbors missing in adata)
    :param iso_threshold: max number of batches per label for label to be considered as
        isolated, if iso_threshold is integer.
        If iso_threshold=None, consider minimum number of batches that labels are present in
    :param return_all: return scores for all isolated labels instead of aggregated mean
    :param verbose:
    :return:
        Mean of scores for each isolated label
        or dictionary of scores for each label if `return_all=True`
    """
    scores = {}
    isolated_labels = get_isolated_labels(
        adata,
        label_key,
        batch_key,
        iso_threshold,
        verbose
    )

    for label in isolated_labels:
        score = score_isolated_label(
            adata,
            label_key,
            label,
            embed,
            cluster,
            verbose=verbose
        )
        scores[label] = score
    scores = pd.Series(scores)

    if return_all:
        return scores

    return scores.mean()


def score_isolated_label(
        adata,
        label_key,
        isolated_label,
        embed,
        cluster=True,
        iso_label_key='iso_label',
        verbose=False
):
    """
    Compute label score for a single label
    :param adata: anndata object
    :param label_key: key in adata.obs of isolated label type (usually cell label)
    :param isolated_label: value of specific isolated label e.g. cell type/identity annotation
    :param embed: embedding to be passed to opt_louvain, if adata.uns['neighbors'] is missing
    :param cluster: if True, compute clustering-based F1 score, otherwise compute
        silhouette score on grouping of isolated label vs all other remaining labels
    :param iso_label_key: name of key to use for cluster assignment for F1 score or
        isolated-vs-rest assignment for silhouette score
    :param verbose:
    :return:
        Isolated label score
    """
    adata_tmp = adata.copy()

    def max_f1(adata, label_key, cluster_key, label, argmax=False):
        """cluster optimizing over largest F1 score of isolated label"""
        obs = adata.obs
        max_cluster = None
        max_f1 = 0
        for cluster in obs[cluster_key].unique():
            y_pred = obs[cluster_key] == cluster
            y_true = obs[label_key] == label
            f1 = f1_score(y_pred, y_true)
            if f1 > max_f1:
                max_f1 = f1
                max_cluster = cluster
        if argmax:
            return max_cluster
        return max_f1

    if cluster:
        # F1-score on clustering
        opt_louvain(
            adata_tmp,
            label_key,
            cluster_key=iso_label_key,
            label=isolated_label,
            use_rep=embed,
            function=max_f1,
            verbose=False,
            inplace=True
        )
        score = max_f1(adata_tmp, label_key, iso_label_key, isolated_label, argmax=False)
    else:
        # AWS score between label
        adata_tmp.obs[iso_label_key] = adata_tmp.obs[label_key] == isolated_label
        score = silhouette(adata_tmp, iso_label_key, embed)

    del adata_tmp

    if verbose:
        print(f"{isolated_label}: {score}")

    return score


def get_isolated_labels(adata, label_key, batch_key, iso_threshold, verbose):
    """
    Get labels that are isolated depending on the number of batches
    :param adata: anndata object
    :param label_key: column in adata.obs
    :param batch_key: column in adata.obs
    :param iso_threshold: Maximum number of batches per label for label to be considered as isolated.
    :param verbose:
    """

    tmp = adata.obs[[label_key, batch_key]].drop_duplicates()
    batch_per_lab = tmp.groupby(label_key).agg({batch_key: "count"})

    # threshold for determining when label is considered isolated
    if iso_threshold is None:
        iso_threshold = batch_per_lab.min().tolist()[0]

    if verbose:
        print(f"isolated labels: no more than {iso_threshold} batches per label")

    labels = batch_per_lab[batch_per_lab[batch_key] <= iso_threshold].index.tolist()
    if len(labels) == 0 and verbose:
        print(f"no isolated labels with less than {iso_threshold} batches")

    return labels


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