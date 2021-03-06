from sklearn.metrics import accuracy_score
from anndata import AnnData
from typing import Dict, Iterable


def accuracy_paired_omics_per_cell_type(
        adata: AnnData,
        bc_list1: Iterable[str],
        bc_list2: Iterable[str],
        omic_layer: str,
        variable: str,
        cell_type: str,
        percent: bool = False
) -> Dict[str, float]:
    """
    will match cell barcode from paired measurement for 2 layers. 
    I will return the dict of ratio of cells for which the RNA and ATAC barcode end up in the same cluster.
    But the ratios are per cell types
    
    Parameters
    ----------
    adata
        coembed multiomic object.
    bc_list1
        RNA matching barcodes.
    bc_list2
        ATAC matching barcodes.
    variable
        cell clustering obs variable.
    omic_layer
        obs variable containing the batch/omic layer of origin.
    cell_type
        obs variable containing the ground truth cell type.
    percent
        If True returns the percentage. If False returns the ratio.
    
    Returns
    -------
    accuracy: dict of float ratio of cells for which the barcodes end up in the same barcodes per cell type
    """

    # extract important informations from the adata.obs
    df = adata.obs[[omic_layer, variable, cell_type]]
    # split RNA and ATAC cells in 2 dataframes
    omic_layer_variable = list(set(df[omic_layer]))
    df_atac = df[df[omic_layer] == omic_layer_variable[0]]
    df_rna = df[df[omic_layer] == omic_layer_variable[1]]
    df_rna = df_rna.loc[bc_list1, :]
    df_atac = df_atac.loc[bc_list2, :]
    cell_type_dict = {}
    for current_cell_type in sorted(set(df[cell_type])):
        df_rna_cell_type = df_rna[df_rna[cell_type] == current_cell_type]
        df_atac_cell_type = df_atac[df_atac[cell_type] == current_cell_type]
        # get the accuracy
        accuracy = accuracy_score(df_rna_cell_type[variable], df_atac_cell_type[variable])
        if percent:
            accuracy = accuracy * 100
        cell_type_dict[current_cell_type] = accuracy

    return cell_type_dict
