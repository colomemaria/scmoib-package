import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from anndata import AnnData
from typing import Union, Tuple


def node_distr(
        adata: AnnData,
        title: Union[str, None] = None,
        figsize: Tuple[float, float] = (10, 5),
        show: bool = True,
        save: Union[str, None] = None
) -> None:
    """
    Histogram plotting the number of nodes between matching barcodes.

    Parameters
    ----------
    adata
        Annotated data matrix with calculated shortest paths in uns part.
    title
        Set the title for the plot.
    figsize
        Width and height in inches.
    show
        Show the plot.
    save
        Save the plot.
    """

    v1, v2 = np.unique(adata.uns['metrics']['nodes_count'], return_counts=True)
    fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(x=v1, y=v2, ax=ax)
    ax.set(xlabel='Number of nodes between barcodes', ylabel='Number of cell pairs', title=title)
    ax.set_xticklabels(v1, rotation=-90)

    if save:
        plt.savefig(save)
    if show:
        plt.show()
    else:
        plt.close()
