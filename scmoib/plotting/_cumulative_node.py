import matplotlib.pyplot as plt
import seaborn as sns


def cumulative_node(adata, title=None, figsize=(10, 10)):
    fig, ax = plt.subplots(figsize=figsize)
    sns.ecdfplot(adata.uns['node_metrics']['nodes_count'])
    sns.kdeplot(adata.uns['node_metrics']['nodes_count'], cumulative=True, color='r', linestyle="dashed")
    ax.set(xlabel='Number of nodes between matching barcodes',
           ylabel='Proportion of matching barcodes within a given number of nodes', 
           title=title)