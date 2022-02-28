import scanpy as sc
#import anndata as ad 
from anndata import AnnData
from typing import Union, Optional, Sequence, Any, Mapping, Tuple, Callable, List, Iterable

def umap_barcodes(
        adata: AnnData,
        bc_list1: Iterable[str],
        bc_list2: Iterable[str],
        basis: str = 'umap',
        color: Union[str, Sequence[str], None] = None,
        edges: bool = True,
        edges_width: float = 0.1
) -> None:
    """

    Parameters
    ----------
    adata
        Annotated data matrix
    bc_list1
        List of barcodes
    bc_list2
        List of barcodes matched with bc_list1
    basis
        Name of the `obsm` basis to use.
    color

    edges
        If True shows the edges between matched barcodes.
    edges_width
        Width of edges.
    """
    bc_list = list(adata.obs.index)
    edge_list = list(zip([bc_list.index(i) for i in bc_list1], [bc_list.index(j) for j in bc_list2]))
    sc.plotting._tools.scatterplots.embedding(adata, basis=basis, color=color, edges=edges, edges_width=edges_width, edge_list=edge_list)

