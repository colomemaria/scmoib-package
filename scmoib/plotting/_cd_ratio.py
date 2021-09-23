import matplotlib.pyplot as plt
from typing import Union
from pandas import DataFrame


def cd_ratio(df: DataFrame, show: bool = True, save: Union[str, None] = None) -> None:
    """
    Barplot for the connected / disconnected barcode ratio for each method.
    Parameters
    ----------
    df
        Metrics dataframe
    show
        Show the plot
    save
        Save the plot
    """
    tdf = df.copy()
    tdf['Connected'] = 1 - tdf['disc_ratio']
    tdf.rename(columns={'disc_ratio': 'Disconnected'}, inplace=True)
    ax = tdf[['Disconnected', 'Connected']].plot(kind='bar', stacked=True)
    ax.legend(bbox_to_anchor=(1.0, 1.0))
    if save:
        plt.savefig(save)
    if show:
        plt.show()
    else:
        plt.close()
