from pandas import DataFrame
from typing import Tuple, Union
from sklearn.preprocessing import MinMaxScaler
import seaborn as sns
import matplotlib.pyplot as plt


def metrics_scatterplot(
        df: DataFrame,
        figsize: Tuple[float, float] = (5, 5),
        title: Union[str, None] = None,
        save: Union[str, None] = None
) -> None:
    """
    Draw a scatter plot for metric values. Missing dot means that metric value is not available.

    Parameters
    ----------
    df
        Dataframe with metric values.
    figsize
        Width and height in inches.
    title
        Title of the plot.
    save
        The name of the file.
    """
    scaler = MinMaxScaler()
    temp_df = DataFrame(scaler.fit_transform(df), columns=df.columns)
    temp_df.index = df.index.copy()
    temp_df = temp_df.stack().reset_index()
    temp_df.rename(columns={'level_0': 'method', 'level_1': 'metric', 0: 'value'}, inplace=True)
    temp_df.method, temp_df.metric = temp_df.method.astype('category'), temp_df.metric.astype('category')
    fig, ax = plt.subplots(figsize=figsize)

    metrics = list(set(temp_df.metric))
    palette = sns.color_palette(n_colors=len(metrics))
    colors_list = []
    for i in temp_df['metric']:
        colors_list.append(palette[metrics.index(i)])
    temp_df['color'] = colors_list

    sns.scatterplot(
        x=temp_df.metric,
        y=temp_df.method,
        size=temp_df.value,
        c=temp_df.color,
        legend=False,
        ax=ax,
    )
    ax.set_title(title)

    if save:
        fig.savefig(save)
