from anndata import AnnData
from scib.metrics.isolated_labels import isolated_labels


def isolated_labels_score(
        adata: AnnData,
        cell_label: str,
        batch_key: str,
        embed: str,
):
    """
    Score how well labels of isolated labels are distiguished in the dataset by either
        1. clustering-based approach F1 score, or
        2. average-width silhouette score (ASW) on isolated label vs all other labels
    Parameters
    ----------
    adata
        Annotated data matrix
    batch_key
        Obs variable containing batch annotation
    cell_label
        Obs variable containing cell type annotation
    embed
        Obsm variable name
    Returns
    -------
    il_score_sil, il_score_clus
    """
    
    il_score_sil = isolated_labels(adata, 
                                   label_key=cell_label, 
                                   batch_key=batch_key, 
                                   embed=embed,
                                   cluster=False, 
                                   verbose=False)
    
    il_score_clus = isolated_labels(adata, 
                                    label_key=cell_label, 
                                    batch_key=batch_key, 
                                    embed=embed,
                                    cluster=True, 
                                    verbose=False)
    
    return il_score_sil, il_score_clus
