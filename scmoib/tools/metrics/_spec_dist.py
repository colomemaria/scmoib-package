from scipy.stats import gaussian_kde
from anndata import AnnData
import numpy as np


def spec_dist(
        adata: AnnData, 
        n_metr: int = 10, 
        norm: bool = True
) -> float:
    """
    Calculate our special distance based on the shortest path statistics.
    This metric is normalized by the ratio of connected barcodes.
    """
    try:
        data = adata.uns['metrics']
        kde = gaussian_kde(data['nodes_count'])
        res = kde.integrate_box_1d(low=-np.inf, high=n_metr)
        if norm:
            return res * data['conn_ratio'] 
        return res
    except KeyError:
        print("Shortest path statistics is not available. \n",
              "Please, check if you ran the nodes_metrics method or used unintegrated data")
        return
