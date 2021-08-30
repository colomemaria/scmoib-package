import matplotlib.pyplot as plt
import seaborn as sns


def cumulative_node(adata, title=None):
    fig, ax = plt.subplots(figsize=(5, 5))
    sns.ecdfplot(adata.uns['nodes_count'])
    sns.kdeplot(adata.uns['nodes_count'], cumulative=True, color='r', linestyle="dashed")
    ax.set(xlabel='Number of nodes between matching barcodes',
           ylabel='Proportion of matching barcodes within a given number of nodes', 
           title=title)