from typing import Union, Iterable, Dict

import pandas as pd
from anndata import AnnData

from . import metrics


class MetricsCalculator:
    """
    The MetricsCalculator object calculates and stores different metrics for the data integration benchmarking.
    """

    def __init__(self):
        self.metrics = {}

    def __check_key(self, key):
        """
        Internal key checking function.
        """
        if key not in self.metrics.keys():
            self.metrics[key] = {}

    def get_df(
            self,
            filter_values: Union[Dict[str, Iterable[str]], None] = None
    ) -> pd.DataFrame:
        """
        Return metrics information as DataFrame.

        Parameters
        ----------
        filter_values
            Dictionary of lists containing non-appropriate metrics for different data integration methods.
            If None returns raw metrics dataframe.
        """
        if filter_values:
            temp_metrics = self.metrics.copy()
            for key in filter_values.keys():
                for val in filter_values[key]:
                    temp_metrics[key][val] = float('nan')
            return pd.DataFrame(temp_metrics).transpose()
        return pd.DataFrame(self.metrics).transpose()

    def node_metrics(
            self,
            adata: AnnData,
            adata_id: str,
            bc_list1: Iterable[str],
            bc_list2: Iterable[str],
            cell_type: str,
            n_jobs: Union[int, None] = None
    ) -> None:
        """
        Calculate node metrics using shortest paths between matching barcodes.
        
        Parameters
        ----------
        adata
            AnnData object.
        adata_id
            Data ID for metrics dataframe.
        bc_list1
            RNA matching barcodes.
        bc_list2
            ATAC matching barcodes.
        cell_type
            obs variable containing the ground truth cell type.
        n_jobs
            The number of jobs to use for the computation. None means 1.
        """
        self.__check_key(adata_id)

        res = metrics.node_metrics(adata, bc_list1, bc_list2, cell_type, n_jobs=n_jobs)
        self.metrics[adata_id]['conn_ratio'] = res[2]

    def spec_dist(
            self,
            adata: AnnData,
            adata_id: str,
            n_metr: int = 10,
            norm: bool = True
    ):
        """
        Calculate our special distance based on the shortest path statistics.
        This metric is normalized by the ratio of connected barcodes.
        """

        self.__check_key(adata_id)
        res = metrics.spec_dist(adata, n_metr=n_metr, norm=norm)
        self.metrics[adata_id][f'spec_dist_{n_metr}'] = res

    def silhouette(
            self,
            adata: AnnData,
            adata_id: str,
            cell_label: str,
            embed: str,
            metric: str = 'euclidean'
    ) -> None:
        """
        Calculate silhouette metrics.
        
        Parameters
        ----------
        adata
            Annotated data matrix.
        adata_id
            Data ID for metrics dataframe.
        cell_label
            obs variable containing cell type information.
        embed
            obsm variable name.
        metric
            Metric type
        """
        self.__check_key(adata_id)
        self.metrics[adata_id]['sil_global'] = metrics.silhouette(adata,
                                                                  cell_label,
                                                                  embed,
                                                                  metric=metric)

    def silhouette_batch(
            self,
            adata: AnnData,
            adata_id: str,
            batch_key: str,
            cell_label: str,
            embed: str,
            metric: str = 'euclidean',
    ) -> None:
        """
            Calculate silhouette metrics.

            Parameters
            ----------
            adata
                Annotated data matrix.
            adata_id
                Data ID for metrics dataframe.
            batch_key
                obs variable containing batch information.
            cell_label
                obs variable containing cell type information.
            embed
                obsm variable name
            metric
                Metric type.
        """
        self.__check_key(adata_id)
        self.metrics[adata_id]['sil_clus'] = metrics.silhouette_batch(adata,
                                                                      batch_key,
                                                                      cell_label,
                                                                      embed,
                                                                      metric=metric)

    def isolated_labels(
            self,
            adata: AnnData,
            adata_id: str,
            batch_key: str,
            group_key: str,
            embed: str,
    ) -> None:
        """
            Calculate silhouette metrics.

            Parameters
            ----------
            adata
                Annotated data matrix.
            adata_id
                Data ID for metrics dataframe.
            batch_key
                obs variable containing batch information.
            group_key
                obs variable containing cell type information.
            embed
                obsm variable name
        """
        self.__check_key(adata_id)
        self.metrics[adata_id]['il_sil'] = metrics.isolated_labels(adata,
                                                                   group_key,
                                                                   batch_key,
                                                                   embed=embed,
                                                                   iso_threshold=1,
                                                                   cluster=False)
        self.metrics[adata_id]['il_clus'] = metrics.isolated_labels(adata,
                                                                    group_key,
                                                                    batch_key,
                                                                    embed=embed,
                                                                    iso_threshold=1,
                                                                    cluster=True)

    def ari(
            self,
            adata: AnnData,
            adata_id: str,
            group1: str,
            group2: str
    ) -> None:
        """
        Calculate Adjusted Rand Index (ARI).
        
        Parameters
        ----------
        adata
            Annotated data matrix.
        adata_id
            Data ID for the metrics dataframe.
        group1
            ground-truth cluster assignments (e.g. cell type labels).
        group2
            "predicted" cluster assignments.
        """
        self.__check_key(adata_id)
        self.metrics[adata_id]['ARI'] = metrics.ari(adata, group1, group2)

    def nmi(
            self,
            adata: AnnData,
            adata_id: str,
            group1: str,
            group2: str,
            method: str = "arithmetic",
            nmi_dir: Union[str, None] = None
    ) -> None:
        """
        Calculate Normalized Mutual Info score (NMI).
        
        Parameters
        ----------
        adata
            Annotated data matrix
        adata_id
            Data ID for metrics dataframe.
        group1
            Ground-truth cluster assignments (e.g. cell type labels).
        group2
            "Predicted" cluster assignments.
        method
            NMI method.
        nmi_dir
            path to the directory with compiled NMI code.
        """
        self.__check_key(adata_id)
        self.metrics[adata_id]['NMI'] = metrics.nmi(adata, group1, group2, method=method, nmi_dir=nmi_dir)

    def ami(
            self,
            adata: AnnData,
            adata_id: str,
            group1: str,
            group2: str
    ) -> None:
        """
        Calculate Adjusted Mutual Info score (AMI).
        
        Parameters
        ----------
        adata:
            Annotated data matrix.
        adata_id
            Data ID for metrics dataframe.
        group1
            Ground-truth cluster assignments (e.g. cell type labels).
        group2
            "Predicted" cluster assignments.
        """
        self.__check_key(adata_id)
        self.metrics[adata_id]['AMI'] = metrics.ami(adata, group1, group2)

    def homogeneity(
            self,
            adata: AnnData,
            adata_id: str,
            group1: str,
            group2: str
    ) -> None:
        """
        Calculate Homogeneity score.
        
        Parameters
        ----------
        adata
            Annotated data matrix.
        adata_id
            Data ID for metrics dataframe.
        group1
            # TODO
        group2
            # TODO
        """
        self.__check_key(adata_id)
        self.metrics[adata_id]['homogeneity'] = metrics.homogeneity(adata, group1, group2)

    def _average_dist_between_matching_barcodes(
            self,
            adata: AnnData,
            adata_id: str,
            bc_list1: Iterable[str],
            bc_list2: Iterable[str],
            metric: str = 'cosine',
            cell_type: str = None,
            absolute: bool = True
    ) -> None:
        """
        Parameters
        ----------
        adata
            coembed multiomic object.
        adata_id
            Data ID for the metrics DataFrame
        bc_list1
            RNA matching barcodes
        bc_list2
            ATAC matching barcodes
        metric
            Metric name ('euclidean' or 'cosine') # TODO
        cell_type
            obs variable containing the ground truth cell type
        absolute
            If True computes the absolute value of distances. Else computes normalized distances.
            """
        result = metrics.average_distance_between_matching_barcodes(adata,
                                                                    bc_list1,
                                                                    bc_list2,
                                                                    metric=metric,
                                                                    cell_type=cell_type,
                                                                    absolute=absolute)

        if isinstance(result, dict):
            adata.uns['average_dist_per_cluster'] = result
        else:
            self.__check_key(adata_id)
            self.metrics[adata_id]['pairwise_distance'] = result

    def _accuracy_paired_omics(
            self,
            adata: AnnData,
            adata_id: str,
            bc_list1: Iterable[str],
            bc_list2: Iterable[str],
            omic_layer: str,
            variable: str,
            percent: bool = False
    ) -> None:
        """
        Will match cell barcode from paired measurement for 2 layers.
        I will return the ratio of cells for which the RNA and ATAC barcode end up in the same cluster.

        Parameters
        ----------
        adata
            coembed multiomic object.
        adata_id
            Data ID for the metrics dataframe.
        bc_list1
            RNA matching barcodes.
        bc_list2
            ATAC matching barcodes.
        variable
            Cell clustering obs variable.
        omic_layer
            obs variable containing the batch/omic layer of origin.
        percent
            If True returns the percentage. If False returns the ratio.

        """

        self.__check_key(adata_id)
        self.metrics[adata_id]['accuracy'] = metrics.accuracy_paired_omics(adata,
                                                                           bc_list1,
                                                                           bc_list2,
                                                                           omic_layer,
                                                                           variable,
                                                                           percent=percent)

    def _accuracy_paired_omics_per_cell_type(
            self,
            adata: AnnData,
            bc_list1: Iterable[str],
            bc_list2: Iterable[str],
            omic_layer: str,
            variable: str,
            cell_type: str,
            percent: bool = False
    ) -> None:

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
        """
        adata.uns['acc_cell_type'] = metrics.accuracy_paired_omics_per_cell_type(adata,
                                                                                 bc_list1,
                                                                                 bc_list2,
                                                                                 omic_layer,
                                                                                 variable,
                                                                                 cell_type,
                                                                                 percent=percent)

    def _graph_connectivity(
            self,
            adata: AnnData,
            adata_id: str,
            label_key: str
    ) -> None:
        """
        Metric that quantifies how connected the subgraph corresponding to each batch cluster is.

        Parameters
        ----------
        adata
            Annotated data matrix.
        adata_id:
            Data ID for the metrics DataFrame.
        label_key
            Batch cluster key.
        """
        self.__check_key(adata_id)
        self.metrics[adata_id][f'graph connectivity {label_key}'] = metrics.graph_connectivity(adata, label_key)
