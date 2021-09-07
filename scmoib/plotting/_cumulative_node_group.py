import matplotlib.pyplot as plt
import seaborn as sns


def cumulative_node_group(adata_list, legend, title=None, figsize=(10, 10)):
    fig, ax = plt.subplots(figsize=figsize)
    for i in adata_list:
        sns.kdeplot(i.uns['node_metrics']['nodes_count'], cumulative=True)
    ax.set(xlabel='Number of nodes between matching barcodes',
           ylabel='Proportion of matching barcodes within a given number of nodes', 
           title=title)
    plt.legend(legend)