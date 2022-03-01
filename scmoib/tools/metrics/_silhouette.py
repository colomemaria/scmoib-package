# This original version of this code was written for the scIB project
# For more information see: https://github.com/theislab/scib
# Paper to cite for this code : https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2
# M. D. Luecken, M. Bu ̈ttner, K. Chaichoompu, A. Danese, M. Interlandi, M. F. Mueller, D. C. Strobl, L. Zappia,
# M. Dugas, M. Colome ́-Tatche ́, and F. J. Theis. Benchmarking atlas-level data integration in single-cell genomics.
# https://www.nature.com/articles/s41592-021-01336-8

import pandas as pd
from anndata import AnnData
from sklearn.metrics.cluster import silhouette_samples, silhouette_score


def silhouette(
        adata: AnnData,
        group_key: str,
        embed: str,
        metric: str = 'euclidean',
        scale: bool = True
) -> float:
    """
    Wrapper for sklearn silhouette function values range from [-1, 1] with
        1 being an ideal fit
        0 indicating overlapping clusters and
        -1 indicating misclassified cells

    Parameters
    ----------
    adata
        Annotated data matrix.
    group_key
        obs variable containing cell type information.
    embed
        obsm variable name.
    metric
        metric type.
    scale
        Scale between 0 (worst) and 1 (best).
    Returns
    -------
    asw
        Silhouette score.
    """
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f'{embed} not in obsm')
    asw = silhouette_score(
        X=adata.obsm[embed],
        labels=adata.obs[group_key],
        metric=metric
    )
    if scale:
        asw = (asw + 1) / 2
    return asw


def silhouette_batch(
        adata: AnnData,
        batch_key: str,
        group_key: str,
        embed: str,
        metric: str = 'euclidean',
        return_all: bool = False,
        scale: bool = True,
        verbose: bool = False
):
    """
    Absolute silhouette score of batch labels subsetted for each group.
    Parameters
    ----------
    adata
        Annotated data matrix.
    batch_key
        obs variable containing batch information.
    group_key
        obs variable containing cell type information.
    embed
        obsm variable name.
    metric
        metric type.
    return_all
        Return all outputs
    scale
        Scale between 0 (worst) and 1 (best).
    verbose
        Display scores per cell type.
    Returns
    -------
    asw
        Mean silhouette score.
    sil_means
        Mean silhouette per group.
    sil_all
        Absolute silhouette scores per group label.
    """
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f'{embed} not in obsm')

    sil_all = pd.DataFrame(columns=['group', 'silhouette_score'])

    for group in adata.obs[group_key].unique():
        adata_group = adata[adata.obs[group_key] == group]
        n_batches = adata_group.obs[batch_key].nunique()

        if (n_batches == 1) or (n_batches == adata_group.shape[0]):
            continue

        sil_per_group = silhouette_samples(
            adata_group.obsm[embed],
            adata_group.obs[batch_key],
            metric=metric
        )

        # take only absolute value
        sil_per_group = [abs(i) for i in sil_per_group]

        if scale:
            # scale s.t. highest number is optimal
            sil_per_group = [1 - i for i in sil_per_group]

        sil_all = pd.concat([sil_all,
                             pd.DataFrame({
                                 'group': [group] * len(sil_per_group),
                                 'silhouette_score': sil_per_group
                             })]
                           )

    sil_all = sil_all.reset_index(drop=True)
    sil_means = sil_all.groupby('group').mean()
    asw = sil_means['silhouette_score'].mean()

    if verbose:
        print(f'mean silhouette per cell: {sil_means}')

    if return_all:
        return asw, sil_means, sil_all

    return asw
