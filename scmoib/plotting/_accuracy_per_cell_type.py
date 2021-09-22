import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Union
from anndata import AnnData


def accuracy_per_cell_type(
        adata: AnnData,
        accuracy: dict,
        cell_type: str,
        display_value: bool = True,
        show: Optional[bool] = True,
        save: Union[str, None] = None
) -> None:
    """
    Accuracy bar plot for each cell type.
    Parameters
    ----------
    adata
        Annotated data matrix
    accuracy
        Dictionary with accuracy score for each cell type
    cell_type
        Obs variable, containing cell type information
    display_value
        Attach a text label above each bar on the plot, displaying its height
    show
        Display plot
    save
        Save plot
    """
    labels = sorted(accuracy.keys())
    key = cell_type + '_colors'
    if key in adata.uns.keys():
        colors = adata.uns[key]
    else:
        colors = sns.color_palette(n_colors=len(labels))
        adata.uns[key] = colors.as_hex()
    values = []
    for label in labels:
        values.append(round(accuracy[label], ndigits=2))
    x = np.arange(len(labels))

    f, ax = plt.subplots(figsize=(18, 9))  # set the size that you'd like (width, height)

    value_plot = plt.bar(x, values, color=colors)
    # Create names on the x-axis
    # ax.legend(fontsize = 14)
    plt.xticks(x, labels, rotation=90)
    ax.set_ylabel('Percentage of matching barcodes clustered together')
    ax.set_title('Percentage of matching barcodes per cell type')
    if display_value:
        # """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in value_plot:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

    if save:
        plt.savefig(save, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close()
