import numpy as np
import matplotlib.pyplot as plt


def accuracy_per_cell_type(adata,
                           accuracy,
                           cell_type,
                           display_value=True,
                           show=True,
                           save=None):
    labels = sorted(accuracy.keys())
    colors = adata.uns[cell_type + '_colors']
    values = []
    for cell_type in labels:
        values.append(round(accuracy[cell_type], ndigits=2))
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
