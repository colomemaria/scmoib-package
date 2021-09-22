import matplotlib.pyplot as plt
import seaborn as sns
from anndata import AnnData
from typing import Sequence, Union, Tuple


def cumulative_node_group(
        adata_list: Sequence[AnnData],
        legend: Tuple[str],
        title: Union[str, None] = None,
        figsize: Tuple[float, float] = (10, 10),
        show: bool = True,
        save: Union[str, None] = None
) -> None:
    """
    Plot the cumulative distribution of the number of connected barcodes for each dataset.
    Parameters
    ----------
    adata_list
        List of annotated data matrices with with calculated shortest path distances in uns parts.
    legend
        Set the legend.
    title
        Set the title for the plot.
    figsize
        Width, height in inches.
    show
        Show the plot.
    save
        Save the plot.
    """
    fig, ax = plt.subplots(figsize=figsize)
    for i in adata_list:
        sns.kdeplot(i.uns['node_metrics']['nodes_count'], cumulative=True)
    ax.set(xlabel='Number of nodes between matching barcodes',
           ylabel='Proportion of matching barcodes within a given number of nodes',
           title=title)
    ax.legend(legend, bbox_to_anchor=(1, 1))
    if save:
        plt.savefig(save)
    if show:
        plt.show()
    else:
        plt.close()
