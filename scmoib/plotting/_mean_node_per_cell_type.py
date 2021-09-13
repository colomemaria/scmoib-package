import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def mean_node_per_cell_type(adata, title=None, figsize=(10, 5)):
    cell_dist = adata.uns['node_metrics']['mean_nodes_per_cell_type']
    fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(x=0, y=1, data=pd.DataFrame(cell_dist.items()))
    ax.set(xlabel='Cell type', ylabel='Mean number of nodes between matching barcodes', title=title)
    ax.set_xticklabels(cell_dist.keys(), rotation=-90)
