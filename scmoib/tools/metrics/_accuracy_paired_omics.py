from sklearn.metrics import accuracy_score


def accuracy_paired_omics(adata, bc_list1, bc_list2, omic_layer, variable, percent=False):
    """
    will match cell barcode from paired measurement for 2 layers. 
    I will return the ratio of cells for which the RNA and ATAC barcode end up in the same cluster.
    
    Parameters
    ----------
    
    adata : coembed multiomic object
    bc_list1: RNA matching barcodes
    bc_list2: ATAC matching barcodes
    variable : cell clustering obs variable
    omic_layer : obs variable containing the batch/omic layer of origin
    percent=True  return percentage. if false, return ratio
    
    Returns
    -------
    
    accuracy: float ratio of cells for which the barcodes end up in the same barcodes
    
    """
    
    # extract important informations from the adata.obs
    df = adata.obs[[omic_layer, variable]]
    omic_layer_variable = list(set(df[omic_layer]))
    df_atac = df[df[omic_layer]==omic_layer_variable[0]]
    df_rna = df[df[omic_layer]==omic_layer_variable[1]]

    df_rna.reindex(index=bc_list1, copy=False)
    df_atac.reindex(index=bc_list2, copy=False)
    
    
    # get the accuracy
    accuracy = accuracy_score(df_rna[variable], df_atac[variable])
    if percent==True:
        accuracy=accuracy*100
    return accuracy