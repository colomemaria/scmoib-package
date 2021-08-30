import numpy as np
import pandas as pd
from .._metrics_stuff import metrics_paired_data as mpd
from .._metrics_stuff import metrics as me
from .._metrics_stuff import dijkstra


class MetricsCalculator:
    def __init__(self):
        self.metrics = {}
        
    def __check_key(self, key):
        if key not in self.metrics.keys():
            self.metrics[key] = {}
        
    def node_metric(self, adata, adata_id, bc_list1, bc_list2, cell_type=None, n_jobs=1):
        self.__check_key(adata_id) 
        results = dijkstra.run_dijkstra(adata, bc_list1, bc_list2, n_jobs)
        tmp_res = list(zip(*results))
        dists = np.array(tmp_res[0])
        num_inf = np.where(dists == float('inf'))[0].shape[0]
        paths = list(np.array(tmp_res[1], dtype=object)[np.where(dists != float('inf'))])
        if num_inf == dists.shape[0]:
            nodes_count = []
            mean_nodes = float('inf')
        else:
            nodes_count = list(map(lambda x: len(x) - 2, paths))
            mean_nodes = np.array(nodes_count).mean()
        self.metrics[adata_id]['num_inf'] = num_inf
        self.metrics[adata_id]['mean_nodes'] = mean_nodes
        self.metrics[adata_id]['cd_ratio'] = 2 * num_inf / adata.shape[0] 
        adata.uns['num_inf'] = num_inf
        adata.uns['mean_nodes'] = mean_nodes
        adata.uns['dists'] = dists
        adata.uns['nodes_count'] = nodes_count
        
        if cell_type:
            cell_type_dist = {}
            for i in paths:
                key = adata.obs.loc[i[0], cell_type]
                cell_type_dist[key] = cell_type_dist.get(key, []) + [len(i) - 2]

            for i in cell_type_dist.keys():
                cell_type_dist[i] = np.mean(cell_dist[i])
            adata.uns['mean_nodes_per_cell_type'] = cell_type_dist
        
    def silhouette(self, adata, adata_id, batch_key='orig.ident', cell_label='paper.cell.type', embed='X_pca'):
        self.__check_key(adata_id)
        sc1, sc2, sc3, sc4 = mpd.run_silhouette_metrics(adata, 
                                                        batch_key=batch_key, 
                                                        cell_label=cell_label, 
                                                        embed=embed)
        
        self.metrics[adata_id]['sil_global'] = sc1
        self.metrics[adata_id]['sil_clus'] = sc2
        self.metrics[adata_id]['il_score_clus'] = sc3
        self.metrics[adata_id]['il_score_sil'] = sc4
    
    def ari(self, adata, adata_id, group1, group2):
        self.__check_key(adata_id)
        self.metrics[adata_id]['ARI'] = me.ari(adata, group1, group2)
        
    def nmi(self, adata, adata_id, group1, group2, method="arithmetic", nmi_dir=None):
        self.__check_key(adata_id)
        self.metrics[adata_id]['NMI'] = me.nmi(adata, group1, group2, method=method, nmi_dir=nmi_dir)
        
    def accuracy_paired_omics(self, adata, adata_id, omic_layer, variable, cell_name=None, percent=False):
        self.__check_key(adata_id)
        self.metrics[adata_id] = mpd.accuracy_paired_omics(adata, 
                                                           omic_layer, 
                                                           variable, 
                                                           cell_name=cell_name, 
                                                           percent=percent)
    
    def accuracy_paired_omics_per_cell_type(self, adata, adata_id, omic_layer, variable, 
                                            cell_type, cell_name=None, percent=False):
        self.__check_key(adata_id)
        adata.uns['acc_cell_type'] = mpd.accuracy_paired_omics_per_cell_type(adata, 
                                                                             omic_layer, 
                                                                             variable, 
                                                                             cell_type, 
                                                                             cell_name=cell_name, 
                                                                             percent=percent)

    def get_df(self):
        return pd.DataFrame(self.metrics).transpose()