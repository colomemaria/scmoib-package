import scib
from anndata import AnnData
from typing import Union


def nmi(
        adata: AnnData,
        group1: str,
        group2: str,
        method: str = "arithmetic",
        nmi_dir: Union[str, None] = None
) -> float:
    """
    Normalized mutual information NMI based on 2 different cluster assignments `group1` and `group2`
    Parameters
    ----------
    adata
        Annotated data matrix
    group1
        column name of `adata.obs` or group assignment
    group2
        column name of `adata.obs` or group assignment
    method
        NMI implementation:
        'max': scikit method with `average_method='max'`
        'min': scikit method with `average_method='min'`
        'geometric': scikit method with `average_method='geometric'`
        'arithmetic': scikit method with `average_method='arithmetic'`
        'Lancichinetti': implementation by A. Lancichinetti 2009 et al.
        'ONMI': implementation by Aaron F. McDaid et al. (https://github.com/aaronmcdaid/Overlapping-NMI) Hurley 2011
    nmi_dir
        Directory of compiled C code if 'Lancichinetti' or 'ONMI' are specified as `method`.
        Compilation should be done as specified in the corresponding README.
    
    Returns
    -------
    Normalized mutual information (NMI) score
    """

    return scib.metrics.nmi(adata, group1, group2, method, nmi_dir)
