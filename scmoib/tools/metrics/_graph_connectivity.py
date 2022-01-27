from scib.metrics.graph_connectivity import graph_connectivity as scib_gr_conn
from anndata import AnnData


def graph_connectivity(adata: AnnData, label_key: str) -> float:
    """"
    Metric that quantifies how connected the subgraph corresponding to each batch cluster is.

    Parameters
    ----------
    adata
        Annotated data matrix.
    label_key
        Batch cluster key.
    Returns
    -------
    Mean graph connectivity score.
    """
    adata.obs[label_key] = adata.obs[label_key].astype('category')
    return scib_gr_conn(adata, label_key)
