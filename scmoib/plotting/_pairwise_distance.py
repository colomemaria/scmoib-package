import scanpy as sc


def pairwise_distance(adata, metric, cell_type, title=None, rotation=90):
    """
    Parameters
    ----------
    adata: AnnData object
    metric: str
        calculated metrics key
    cell_type: str    
    """
    key = f'{metric}_pairwise_distance_between_matching_barcodes'
    if key not in adata.obs.keys():
        raise KeyError(f'Please compute the {metric} metric for your dataset')
        
    sc.pl.violin(pbmc, keys=key, groupby=cell_type, title=title, rotation=rotation)
