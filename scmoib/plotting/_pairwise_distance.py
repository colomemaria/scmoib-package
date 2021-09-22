import scanpy as sc
from anndata import AnnData
from typing import Union, Optional


def pairwise_distance(
        adata: AnnData,
        cell_type: str,
        metric: str = 'euclidean',
        title: Optional[Union[str, None]] = None,
        rotation: Optional[float] = 90,
        save: Union[str, None] = None
) -> None:
    """
    Violin plot for pairwise distances for each cell type.

    Parameters
    ----------
    adata
        Annotated data matrix with calculated pairwise distances.
    cell_type
        Obs variable containing the cell type annotation.
    metric

    title
        Set the title fir the plot.
    rotation
        Rotation of xtick labels.
    save
        Save the plot
    """
    key = f'{metric}_pairwise_distance_between_matching_barcodes'
    if key not in adata.obs.keys():
        raise KeyError(f'Please compute the {metric} metric for your dataset')
        
    sc.pl.violin(adata, keys=key, groupby=cell_type, title=title, rotation=rotation, save=save)
