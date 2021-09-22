import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from anndata import AnnData
from typing import Union, Tuple


def mean_node_per_cell_type(
        adata: AnnData,
        cell_type: str,
        title: Union[str, None] = None,
        figsize: Tuple[float, float] = (10, 5),
        show: bool = True,
        save: Union[str, None] = None
) -> None:
    """
    Mean node distance barplot for each cell type.
    Parameters
    ----------
    adata
        Annotated data matrix with calculated shortest paths in uns part.
    cell_type
        Obs variable containing cell type annotation.
    title
        Set the title for the plot.
    figsize
        Width and height in inches.
    show
        Show the plot.
    save
        Save the plot.
    """
    key = cell_type + '_colors'
    cell_dist = adata.uns['node_metrics']['mean_nodes_per_cell_type']
    cell_dist = dict(sorted(cell_dist.items(), key=lambda x: x[0]))
    fig, ax = plt.subplots(figsize=figsize)
    if key in adata.uns.keys():
        colors = adata.uns[key]
    else:
        colors = sns.color_palette(n_colors=len(cell_dist.keys()))
        adata.uns[key] = colors.as_hex()
    sns.barplot(x=0, y=1, data=pd.DataFrame(cell_dist.items()), palette=colors)
    ax.set(xlabel='Cell type', ylabel='Mean number of nodes between matching barcodes', title=title)
    ax.set_xticklabels(cell_dist.keys(), rotation=-90)

    if save:
        plt.savefig(save)
    if show:
        plt.show()
    else:
        plt.close()
