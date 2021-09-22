import episcanpy as epi
from anndata import AnnData
from typing import Union


# TODO remove extra parameters
def silhoutte(
        adata: AnnData,
        cluster_annot: str,
        value: str = 'X_pca',
        metric: str = 'euclidean',
        save: Union[str, None] = None,
        palette: object = None,
        key: object = None,
        xlabel: object = None,
        ylabel: object = None,
        title: Union[str, None] = None,
        size: str = 'large',
        alternative_plot: bool = False,
        name_cluster: bool = True,
        name_cluster_pos: str = 'left'
) -> None:
    """
    Parameters
    ----------
    adata
        Annotated data matrix.
    cluster_annot
        obs variable containing cluster annotation.
    value
        Obsm variable.
    metric
        Metric type.
    save
        Save the plot.
    palette
    key
    xlabel
    ylabel
    title
    size
    alternative_plot
    name_cluster
    name_cluster_pos
    """
    epi.tl.silhouette(adata,
                      cluster_annot,
                      value=value,
                      metric=metric)
    epi.pl.silhouette(adata,
                      cluster_annot,
                      value=value,
                      metric=metric,
                      key=key,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      size=size,
                      alternative_plot=alternative_plot,
                      name_cluster=name_cluster,
                      name_cluster_pos=name_cluster_pos,
                      palette=palette,
                      save=save,
                      )
