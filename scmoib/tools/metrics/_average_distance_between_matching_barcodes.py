from sklearn.metrics import pairwise_distances
import numpy as np
from typing import Union, Dict, Iterable
from anndata import AnnData


def __get_pairs(adata, bc_list1, bc_list2):
    bc_list = list(adata.obs.index)
    pairs = list(zip([bc_list.index(i) for i in bc_list1], [bc_list.index(j) for j in bc_list2]))
    return pairs


def __distance_between_matching_barcodes(adata, bc_list1, bc_list2, metric='euclidean', absolute=True):
    distances = pairwise_distances(adata.X, metric=metric)
    distances_mean = distances.mean(axis=0)
    pairs = __get_pairs(adata, bc_list1, bc_list2)
    list_distance_barcodes = [np.nan] * len(adata.obs_names)

    if absolute:
        for position in pairs:
            list_distance_barcodes[position[0]] = distances[position[0], position[1]]
            list_distance_barcodes[position[1]] = distances[position[0], position[1]]
    else:
        for position in pairs:
            list_distance_barcodes[position[0]] = distances[position[0], position[1]] / distances_mean[position[0]]
            list_distance_barcodes[position[1]] = distances[position[0], position[1]] / distances_mean[position[0]]
    adata.obs[f'{metric}_pairwise_distance_between_matching_barcodes'] = list_distance_barcodes


def average_distance_between_matching_barcodes(
        adata: AnnData,
        bc_list1: Iterable[str],
        bc_list2: Iterable[str],
        metric: str = 'euclidean',
        cell_type: Union[str, None] = None,
        absolute: bool = True
) -> Union[float, Dict[str, float]]:
    """
    Computes the global average distance or the average distances for each cell type.

    Parameters
    ----------
    adata
        Annotated data matrix
    bc_list1
        RNA matching barcodes
    bc_list2
        ATAC matching barcodes
    metric
        Metric name ('euclidean' or 'cosine') # TODO
    cell_type
        Obs variable containing cell type annotation
    absolute
        If True computes the absolute value of distances.
        Else computes normalized distances.

    Returns
    -------
    If cell_type is None returns the average distance between all matching barcodes.
    Else returns the dictionary of average distances for each cell type
    """
    __distance_between_matching_barcodes(adata, bc_list1, bc_list2, metric=metric, absolute=absolute)

    if not cell_type:
        average_metric = np.mean(adata.obs[f'{metric}_pairwise_distance_between_matching_barcodes'])
        return average_metric
    else:
        average_dist_per_cluster = {}
        for elem in list(set(adata.obs[cell_type])):
            average_dist_per_cluster[elem] = \
                np.mean(
                    adata[adata.obs[cell_type] == elem, :].obs[f'{metric}_pairwise_distance_between_matching_barcodes'])
        return average_dist_per_cluster
