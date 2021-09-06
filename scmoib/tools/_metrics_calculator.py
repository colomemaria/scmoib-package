import numpy as np
import pandas as pd
import sklearn
import episcanpy as epi
from . import metrics


class MetricsCalculator:
    """
    The MetricsCalculator object calculates and stores different metrics for
    data integration benchmarking.
    """
    def __init__(self):
        self.metrics = {}
        
    def __check_key(self, key):
        if key not in self.metrics.keys():
            self.metrics[key] = {}
            
    def get_df(self):
        """
        Return metrics information as DataFrame.
        """
        return pd.DataFrame(self.metrics).transpose()
        
    def node_metrics(self, adata, adata_id, bc_list1, bc_list2, cell_type, n_jobs=None):
        """
        Calculate node metrics using shortest paths between matching barcodes.
        
        Parameters
        ----------
        adata: AnnData object
            AnnData object
        
        adata_id: str
            Data ID for metrics dataframe.
            
        bc_list1: list
            
        bc_list2: list
        
        cell_type: str
            obs variable containing the ground truth cell type
        
        n_jobs: None or int, optional
            Number of threads for calculating shortest paths
        """
        self.__check_key(adata_id) 
        
        res = metrics.node_metrics(adata, bc_list1, bc_list2, cell_type, n_jobs=n_jobs)
        self.metrics[adata_id]['num_inf'] = res[0]
        self.metrics[adata_id]['mean_nodes'] = res[1]
        
    def silhouette(self, adata, adata_id, batch_key='orig.ident', cell_label='paper.cell.type', embed='X_pca'):
        """
        Calculate silhouette metrics.
        
        Parameters
        ----------
        adata: AnnData object
            AnnData object
        
        adata_id: str
            Data ID for metrics dataframe.
        
        ...    
        """
        self.__check_key(adata_id)
        sc1, sc2, sc3, sc4 = metrics.silhouette(adata, 
                                                batch_key=batch_key, 
                                                cell_label=cell_label, 
                                                embed=embed)
        
        self.metrics[adata_id]['sil_global'] = sc1
        self.metrics[adata_id]['sil_clus'] = sc2
        self.metrics[adata_id]['il_score_clus'] = sc3
        self.metrics[adata_id]['il_score_sil'] = sc4
    
    def ari(self, adata, adata_id, group1, group2):
        """
        Calculate Adjusted Rand Index (ARI)
        
        Parameters
        ----------
        adata: AnnData object
            AnnData object
        
        adata_id: str
            Data ID for metrics dataframe.
        
        group1: str
            ground-truth cluster assignments (e.g. cell type labels)
        group2: str
            "predicted" cluster assignments
        """
        self.__check_key(adata_id)
        self.metrics[adata_id]['ARI'] = metrics.ari(adata, group1, group2)
        
    def nmi(self, adata, adata_id, group1, group2, method="arithmetic", nmi_dir=None):
        """
        Calculate Normalized Mutual Info score (NMI)
        
        Parameters
        ----------
        adata: AnnData object
            AnnData object
        
        adata_id: str
            Data ID for metrics dataframe.
        
        group1: str
            ground-truth cluster assignments (e.g. cell type labels)
        group2: str
            "predicted" cluster assignments
        """
        self.__check_key(adata_id)
        self.metrics[adata_id]['NMI'] = metrics.nmi(adata, group1, group2, method=method, nmi_dir=nmi_dir)
    
    def ami(self, adata, adata_id, group1, group2):
        """
        Calculate Adjusted Mutual Info score (NMI)
        
        Parameters
        ----------
        adata: AnnData object
            AnnData object
        
        adata_id: str
            Data ID for metrics dataframe.
        
        group1: str
            ground-truth cluster assignments (e.g. cell type labels)
        group2: str
            "predicted" cluster assignments
        """
        self.__check_key(adata_id)
        self.metrics[adata_id]['AMI'] = metrics.ami(adata, group1, group2)
        
    def homogeneity(self, adata, adata_id, group1, group2):
        """
        Calculate Homogeneity score (NMI)
        
        Parameters
        ----------
        adata: AnnData object
            AnnData object
        
        adata_id: str
            Data ID for metrics dataframe.
        
        group1: str
        
        group2: str
        """
        self.__check_key(adata_id)
        self.metrics[adata_id]['homogeneity'] = metrics.homogeneity(adata, group1, group2)
    
    def _average_dist_between_matching_barcodes(self, adata, adata_id, cell_names=None, cell_type=None, absolute=True):
        result = metrics.average_distance_between_matching_barcodes(adata, 
                                                                cell_names=cell_names, 
                                                                cell_type=cell_type, 
                                                                absolute=absolute)
        if isinstance(result, float):
            self.__check_key(adata_id)
            self.metrics[adata_id]['pairwise_distance'] = result
        elif isinstance(result, dict):
            adata.uns['average_dist_per_cluster'] = result

    def _accuracy_paired_omics(self, adata, adata_id, omic_layer, variable, cell_name=None, percent=False):
        """
        will match cell barcode from paired measurement for 2 layers. 
        I will return the ratio of cells for which the RNA and ATAC barcode end up in the same cluster.

        Parameters
        ----------
        adata : coembed multiomic object
        variable : cell clustering obs variable
        cell_name : obs variable containing the matching barcodes for the 2 omic layers
        omic_layer : obs variable containing the batch/omic layer of origin
        percent=True  return percentage. if false, return ratio

        Returns
        -------
        accuracy: float ratio of cells for which the barcodes end up in the same barcodes

        """

        self.__check_key(adata_id)
        self.metrics[adata_id] = metrics.accuracy_paired_omics(adata, 
                                                           omic_layer, 
                                                           variable, 
                                                           cell_name=cell_name, 
                                                           percent=percent)
    
    def _accuracy_paired_omics_per_cell_type(self, adata, adata_id, omic_layer, variable, 
                                             cell_type, cell_name=None, percent=False):
        
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
        adata.uns['acc_cell_type'] = metrics.accuracy_paired_omics_per_cell_type(adata, 
                                                                                 omic_layer, 
                                                                                 variable, 
                                                                                 cell_type, 
                                                                                 cell_name=cell_name, 
                                                                                 percent=percent)