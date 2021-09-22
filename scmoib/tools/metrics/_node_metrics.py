from .utils import dijkstra
import numpy as np
from anndata import AnnData
from typing import List, Union, Tuple


def node_metrics(
        adata: AnnData,
        bc_list1: List[str],
        bc_list2: List[str],
        cell_type: str,
        n_jobs: Union[int, None] = None
) -> Tuple[float, float, float]:
    """
    Calculates node metrics using shortest paths between matching barcodes.
        
    Parameters
    ----------
    adata
        AnnData object.
    bc_list1
        RNA matching barcodes.
    bc_list2
        ATAC matching barcodes.
    cell_type
        obs variable containing the ground truth cell type.
    n_jobs
        The number of jobs to use for the computation. None means 1.
    
    Returns
    -------
    num_inf: int
        Number of disconnected barcodes.
    mean_nodes
        Mean node distance between barcodes.
    disc_ratio
        The fraction of disconnected barcodes from the total.
    """
    results = dijkstra.run_dijkstra(adata, bc_list1, bc_list2, n_jobs=n_jobs)
    tmp_res = list(zip(*results))
    dists = np.array(tmp_res[0])
    num_inf = 2 * np.where(dists == float('inf'))[0].shape[0]
    disc_ratio = num_inf / len(adata.obs.index)
    paths = list(np.array(tmp_res[1], dtype=object)[np.where(dists != float('inf'))])

    if num_inf == dists.shape[0]:
        nodes_count = []
        mean_nodes = float('inf')
    else:
        nodes_count = list(map(lambda x: len(x) - 2, paths))
        mean_nodes = np.array(nodes_count).mean()

    cell_type_dist = {}
    for i in paths:
        key = adata.obs.loc[i[0], cell_type]
        cell_type_dist[key] = cell_type_dist.get(key, []) + [len(i) - 2]

    for i in cell_type_dist.keys():
        cell_type_dist[i] = np.mean(cell_type_dist[i])

    metrics = {'num_inf': num_inf, 'mean_nodes': mean_nodes, 'dists': dists, 'nodes_count': nodes_count,
               'mean_nodes_per_cell_type': cell_type_dist, 'disc_ratio': disc_ratio}
    adata.uns['metrics'] = metrics

    return num_inf, mean_nodes, disc_ratio
