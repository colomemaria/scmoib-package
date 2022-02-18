import episcanpy as epi
from anndata import AnnData


def check_anndata(adata: AnnData,
                  cell_label: str,
                  mode: str,
                  use_rep: str = None):
    """
    Prepares data for further analysis depending on used data integration method.
    Generates louvain cluster labels.
    #TODO check the cell bales and barcodes are present
    Parameters
    ----------
    adata
        Anndata object with single cell data
    cell_label
        obs variable containing the ground truth cell type.
    mode
        Available modes: ['raw', 'count', 'comps', 'graph']
    use_rep
        Embedding used for neighbors graph recalculation.
        (only for methods that change components)
    """
    print('Running pre-flight check')
    if mode not in ['raw', 'count', 'comps', 'graph']:
        raise Exception('Functions works only in 4 modes: raw, count, comps, graph')

    if mode == 'raw' or mode == 'count':
        epi.pp.lazy(adata)

    if mode == 'comps':
        epi.pp.neighbors(adata, use_rep=use_rep)
        epi.tl.umap(adata)

    if mode == 'graph':
        epi.tl.umap(adata)
    print('Clustering...')
    epi.tl.getNClusters(adata, adata.obs[cell_label].unique().shape[0])