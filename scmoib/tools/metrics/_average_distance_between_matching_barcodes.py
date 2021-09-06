from sklearn.metrics import pairwise_distances
import networkx as nx
import numpy as np


def __get_pairs(adata, cell_names=None):
    # extract graph tuples for paired cells
    df = adata.obs.copy()
    df

    g = nx.Graph(adata.obsp['connectivities'])
    node_list = []
    for node in g.nodes:
        node_list.append(node)
    df['nodes'] = node_list

    #cell_names = [x.split('_')[0] for x in df.index.tolist()]
    if cell_names == None:
        df['org_cell_name'] = cell_names
    
    df['org_cell_name'] = df[cell_names]
    dico_tuple = {}
    tuples = []
    index = 0
    for name in df['org_cell_name']:
        if name not in dico_tuple.keys():
            dico_tuple[name] = [df['nodes'][index]]
        else:
            dico_tuple[name].append(df['nodes'][index])
            tuples.append(dico_tuple[name])
        index += 1

    return (tuples)


def __distance_between_matching_barcodes(adata, cell_names, absolute=True):
    distances = pairwise_distances(adata.X)
    distances_mean = distances.mean(axis=0)
    pairs = __get_pairs(adata, cell_names=cell_names)
    list_distance_barcodes = [np.nan]*len(adata.obs_names)

    if absolute == True:
        for position in pairs:
            list_distance_barcodes[position[0]] = distances[position[0], position[1]]
            list_distance_barcodes[position[1]] = distances[position[0], position[1]]
    else:
        for position in pairs:
            list_distance_barcodes[position[0]] = distances[position[0], position[1]]/distances_mean[position[0]]
            list_distance_barcodes[position[1]] = distances[position[0], position[1]]/distances_mean[position[0]]
            #return(list_distance_barcodes)
    adata.obs['euclidean_pairwise_distance_between_matching_barcodes'] = list_distance_barcodes

    
def average_distance_between_matching_barcodes(adata, cell_names=None, cell_type=None, absolute=True):
    __distance_between_matching_barcodes(adata, cell_names, absolute)
    
    if cell_type == None:
        average_metric = np.mean(adata.obs['euclidean_pairwise_distance_between_matching_barcodes'])
        return(average_metric)
    else:
        average_dist_per_cluster = {}
        for elem in list(set(adata.obs[cell_type])):
            average_dist_per_cluster[elem] = \
            np.mean(adata[adata.obs[cell_type] == elem,:].obs['euclidean_pairwise_distance_between_matching_barcodes'])
        return average_dist_per_cluster