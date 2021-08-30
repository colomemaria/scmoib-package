import matplotlib.pyplot as plt
import seaborn as sns


def cumulative_node_group(adata_list, legend, title=None):
    fig, ax = plt.subplots(figsize=(5, 5))
    for i in adata_list:
        sns.kdeplot(i.uns['nodes_count'], cumulative=True)
    ax.set(xlabel='Number of nodes between matching barcodes',
           ylabel='Proportion of matching barcodes within a given number of nodes', 
           title=title)
    plt.legend(legend)