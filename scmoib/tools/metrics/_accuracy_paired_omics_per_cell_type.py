def accuracy_paired_omics_per_cell_type(adata, omic_layer, variable, cell_type, cell_name=None, percent=False):
    """
    will match cell barcode from paired measurement for 2 layers. 
    I will return the dict of ratio of cells for which the RNA and ATAC barcode end up in the same cluster.
    But the ratios are per cell types
    
    Parameters
    ----------
    
    adata : coembed multiomic object
    variable : cell clustering obs variable
    cell_name : obs variable containing the matching barcodes for the 2 omic layers
    omic_layer : obs variable containing the batch/omic layer of origin
    cell_type : obs variable containing the ground truth cell type
    percent=True  return percentage. if false, return ratio
    
    Returns
    -------
    
    accuracy: dict of float ratio of cells for which the barcodes end up in the same barcodes per cell type
    
    """
    
    # extract important informations from the adata.obs
    df = adata.obs[[omic_layer, variable, cell_type]]
    if cell_name != None:
        df.index = adata.obs[cell_name]
    else:
        df.index = adata.index

    cell_type_dict = {}
    for current_cell_type in sorted(set(df[cell_type])):
        df_cell_type = df[df[cell_type]==current_cell_type]
        
        
        # split RNA and ATAC cells in 2 dataframes
        omic_layer_variable = list(set(df[omic_layer]))
        df_atac = df_cell_type[df_cell_type[omic_layer]==omic_layer_variable[0]]
        df_rna = df_cell_type[df_cell_type[omic_layer]==omic_layer_variable[1]]
        
        # only keep cells that are present in the 2 tables
        atac_cells = []
        for name in df_atac.index.tolist():
            if name in df_rna.index:
                atac_cells.append(True)
            else:
                atac_cells.append(False)
        
        rna_cells = []
        for name in df_rna.index.tolist():
            if name in df_atac.index:
                rna_cells.append(True)
            else:
                rna_cells.append(False)
        
        df_rna = df_rna[rna_cells].sort_index()
        df_atac = df_atac[atac_cells].sort_index()
    
        # get the accuracy
        x =(df_rna[variable] == df_atac[variable]).tolist()
        Correct_values = x.count(True)
        False_values = x.count(False)
        accuracy = Correct_values/(Correct_values+False_values)
        if percent==True:
            accuracy=accuracy*100
        cell_type_dict[current_cell_type] = accuracy
        
    return cell_type_dict