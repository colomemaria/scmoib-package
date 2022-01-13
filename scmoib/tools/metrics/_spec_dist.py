from scipy.special import ndtr
from scipy.stats import gaussian_kde
from anndata import AnnData

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
        cdf = tuple(ndtr(np.ravel(item - kde.dataset) / kde.factor).mean() for item in np.unique(data))
        if norm:
            return cdf[min(n_metr, len(cdf) - 1)] * data['conn_ratio'] 
        return cdf[min(n_metr, len(cdf) - 1)]
    except KeyError:
        print("Shortest path statistics is not available. \n",
              "Please, check if you ran the nodes_metrics method or used unintegrated data")
        return