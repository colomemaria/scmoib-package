import matplotlib.pyplot as plt
import seaborn as sns
from anndata import AnnData
from typing import Union, Tuple


def cumulative_node(
        adata: AnnData,
        title: Union[str, None] = None,
        figsize: Tuple[float, float] = (10, 10),
        show: bool = True,
        save: Union[str, None] = None
) -> None:
    """

    Parameters
    ----------
    adata
        Annotated data matrix with calculated shortest path distances in uns part.
    title
        Set title for the plot.
    figsize
        Width, height in inches.
    show
        Show the plot.
    save
        Save the plot.
    """
    fig, ax = plt.subplots(figsize=figsize)
    sns.ecdfplot(adata.uns['metrics']['nodes_count'])
    sns.kdeplot(adata.uns['metrics']['nodes_count'], cumulative=True, color='r', linestyle="dashed")
    ax.set(xlabel='Number of nodes between matching barcodes',
           ylabel='Proportion of matching barcodes within a given number of nodes',
           title=title)
    if save:
        plt.savefig(save)
    if show:
        plt.show()
    else:
        plt.close()
