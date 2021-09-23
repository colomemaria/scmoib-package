import matplotlib.pyplot as plt
import seaborn as sns
from anndata import AnnData
from typing import Union, Tuple


def dists_distr(
        adata: AnnData,
        title: Union[str, None] = None,
        bins: Union[str, int] = 'auto',
        figsize: Tuple[float, float] = (10, 5),
        show: bool = True,
        save: Union[str, None] = None
) -> None:
    """
    Historgram plotting the weighted distance between matching barcodes. 
    Parameters
    ----------
    adata
        Annotated data matrix with calculated shortest paths in uns part.
    title
        Set the title for the plot.
    bins
        Generic bin parameter that can be the name of a reference rule, the number of bins, or the breaks of the bins. 
        Passed to numpy.histogram_bin_edges() (https://seaborn.pydata.org/generated/seaborn.histplot.html).
    figsize
        Width and height in inches.
    show
        Show the plot.
    save
        Save the plot.
    """
    fig, ax = plt.subplots(figsize=figsize)
    sns.histplot(data=adata.uns['metrics']['dists'], ax=ax, bins=bins)
    ax.set(xlabel='Weighted distance between barcodes', title=title)
    if save:
        plt.savefig(save)
    if show:
        plt.show()
    else:
        plt.close()
