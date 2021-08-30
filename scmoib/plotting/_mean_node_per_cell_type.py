import matplotlib.pyplot as plt
import seaborn as sns


def mean_node_per_cell_type(adata, title=None):
    cell_dist = adata.uns['mean_nodes_per_cell_type']
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.barplot(x=0, y=1, data=pd.DataFrame(cell_dist.items()))
    ax.set(xlabel='Cell type', ylabel='Mean number of nodes between matching barcodes', title=title)
    ax.set_xticklabels(cell_dist.keys(), rotation=-90);