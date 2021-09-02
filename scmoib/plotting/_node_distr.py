import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def node_distr(adata, title=None):
	"""
	Historgram plotting the number of nodes between matching barcodes. 
	
	"""

    v1, v2 = np.unique(adata.uns['nodes_count'], return_counts=True)
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.barplot(x=v1, y=v2, ax=ax)
    ax.set(xlabel='Number of nodes between barcodes', ylabel='Number of cell pairs', title=title)
    ax.set_xticklabels(v1, rotation=-90);