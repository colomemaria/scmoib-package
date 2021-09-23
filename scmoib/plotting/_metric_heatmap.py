from sklearn.preprocessing import MinMaxScaler
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def metric_heatmap(dataframe, cmap='viridis', scale=False, display_values=True, show=True, save=None):
    """
    for raw we use 'viridis'
    for scaled we use 'plasma'
    """
    if scale:
        # scaling
        scaler = MinMaxScaler()
        index_save = dataframe.index.copy()
        dataframe = pd.DataFrame(scaler.fit_transform(dataframe), columns=dataframe.columns)
        dataframe.index = index_save.copy()

    # plot
    sns.set_style('ticks')
    fig, ax = plt.subplots()
    # the size of A4 paper
    fig.set_size_inches(11.7, 8.27)
    sns.heatmap(dataframe.transpose(copy=True),
                cmap=cmap,
                linewidths=0.2,
                annot=display_values,
                xticklabels=False,
                yticklabels=False,
                square=True)
    plt.yticks(np.arange(0.5, len(dataframe.columns), 1), dataframe.columns)
    plt.xticks(np.arange(0.5, len(dataframe.index), 1), dataframe.index, rotation=90)

    if save:
        fig.savefig(save, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()
