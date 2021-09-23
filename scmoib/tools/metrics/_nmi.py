# This original version of this code was written for the scIB project
# For more information see: https://github.com/theislab/scib
# Paper to cite for this code : https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2
# M. D. Luecken, M. Bu ̈ttner, K. Chaichoompu, A. Danese, M. Interlandi, M. F. Mueller, D. C. Strobl, L. Zappia,
# M. Dugas, M. Colome ́-Tatche ́, and F. J. Theis. Benchmarking atlas-level data integration in single-cell genomics.
# bioRxiv, page 2020.05.22.111161, May 2020. doi: 10.1101/2020.05.22.111161.

from .utils.utils import checkAdata, checkBatch
import sklearn
import pandas as pd
import subprocess
import tempfile
from os import remove
from anndata import AnnData
from typing import Union

ERROR_MESSAGE = 'Please provide the directory of the compiled C code from' \
                'https://sites.google.com/site/andrealancichinetti/mutual3.tar.gz'


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
    nmi_value
        normalized mutual information (NMI)
    """

    checkAdata(adata)

    if isinstance(group1, str):
        checkBatch(group1, adata.obs)
        group1 = adata.obs[group1].tolist()
    elif isinstance(group1, pd.Series):
        group1 = group1.tolist()

    if isinstance(group2, str):
        checkBatch(group2, adata.obs)
        group2 = adata.obs[group2].tolist()
    elif isinstance(group2, pd.Series):
        group2 = group2.tolist()

    if len(group1) != len(group2):
        raise ValueError(f'different lengths in group1 ({len(group1)}) and group2 ({len(group2)})')

    # choose method
    if method in ['max', 'min', 'geometric', 'arithmetic']:
        nmi_value = sklearn.metrics.normalized_mutual_info_score(group1, group2, average_method=method)
    elif method == "Lancichinetti":
        nmi_value = nmi_Lanc(group1, group2, nmi_dir=nmi_dir)
    elif method == "ONMI":
        nmi_value = __onmi(group1, group2, nmi_dir=nmi_dir)
    else:
        raise ValueError(f"Method {method} not valid")

    return nmi_value


def __onmi(group1, group2, nmi_dir=None, verbose=True):
    """
    Based on implementation https://github.com/aaronmcdaid/Overlapping-NMI
    publication: Aaron F. McDaid, Derek Greene, Neil Hurley 2011
    params:
        nmi_dir: directory of compiled C code
    """

    if nmi_dir is None:
        raise FileNotFoundError(ERROR_MESSAGE)

    group1_file = __write_tmp_labels(group1, to_int=False)
    group2_file = __write_tmp_labels(group2, to_int=False)

    nmi_call = subprocess.Popen(
        [nmi_dir + "onmi", group1_file, group2_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT)

    stdout, stderr = nmi_call.communicate()
    if stderr:
        print(stderr)

    nmi_out = stdout.decode()
    if verbose:
        print(nmi_out)

    nmi_split = [x.strip().split('\t') for x in nmi_out.split('\n')]
    nmi_max = float(nmi_split[0][1])

    # remove temporary files
    remove(group1_file)
    remove(group2_file)

    return nmi_max


def nmi_Lanc(group1, group2, nmi_dir="external/mutual3/", verbose=True):
    """
    paper by A. Lancichinetti 2009
    https://sites.google.com/site/andrealancichinetti/mutual
    recommended by Malte
    """

    if nmi_dir is None:
        raise FileNotFoundError(ERROR_MESSAGE)

    group1_file = __write_tmp_labels(group1, to_int=False)
    group2_file = __write_tmp_labels(group2, to_int=False)

    nmi_call = subprocess.Popen(
        [nmi_dir + "mutual", group1_file, group2_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT)

    stdout, stderr = nmi_call.communicate()
    if stderr:
        print(stderr)
    nmi_out = stdout.decode().strip()

    return float(nmi_out.split('\t')[1])


def __write_tmp_labels(group_assignments, to_int=False, delim='\n'):
    """
    write the values of a specific obs column into a temporary file in text format
    needed for external C NMI implementations (onmi and nmi_Lanc functions), because they require files as input
    params:
        to_int: rename the unique column entries by integers in range(1,len(group_assignments)+1)
    """
    if to_int:
        label_map = {}
        i = 1
        for label in set(group_assignments):
            label_map[label] = i
            i += 1
        labels = delim.join([str(label_map[name]) for name in group_assignments])
    else:
        labels = delim.join([str(name) for name in group_assignments])

    clusters = {label: [] for label in set(group_assignments)}
    for i, label in enumerate(group_assignments):
        clusters[label].append(str(i))

    output = '\n'.join([' '.join(c) for c in clusters.values()])
    output = str.encode(output)

    # write to file
    with tempfile.NamedTemporaryFile(delete=False) as f:
        f.write(output)
        filename = f.name

    return filename
