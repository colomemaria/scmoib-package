import pandas as pd
import numpy as np
from plotly.offline import iplot
import plotly.io as pio
from scanpy.plotting._tools.scatterplots import _get_palette


def river_plot(adata,
               source,
               target,
               cell_number=True,
               title='River plot (Sankey Diagram)',
               save=None,
               scale=1):
    """
    cell_number : Bool If True prints the number of cells in each categorie.
    Else, doesn't display it
    scale : float number. above 1 it increase the resolution. below 1 it reudce the resolution
    only matter when saving the plot.
    """
    adata.obs[source] = adata.obs[source].astype('str').astype('category')
    adata.obs[target] = adata.obs[target].astype('str').astype('category')
    df_nodes, df_links = __tool_sankey(adata, source, target, cell_number=cell_number)
    __plot_sankey(df_nodes, df_links,
                  title=title,
                  save=save,
                  scale=scale)


def __tool_sankey(adata, source, target, cell_number=True):
    # extract key_infos in adata
    key_infos = pd.crosstab(adata.obs[target], adata.obs[source])

    ###### NODES ######
    # transform key_infos into the nodes df
    nodes = [['ID', 'Label', 'Color']]
    if cell_number == False:
        label_list = key_infos.columns.tolist() + key_infos.index.tolist()
    else:
        target_cell_nb = pd.crosstab(adata.obs[target], adata.obs[target], margins=True)
        source_cell_nb = pd.crosstab(adata.obs[source], adata.obs[source], margins=True)

        source_names = []
        for n in range(0, len(key_infos.columns.tolist())):
            source_names.append(" n=".join([str(key_infos.columns.tolist()[n]), str(source_cell_nb['All'][n])]))

        target_names = []
        index = 0
        for target_name in key_infos.index.tolist():
            # print(target_name, target_cell_nb['All'][index])
            target_names.append(" n=".join([target_name, str(target_cell_nb['All'][index])]))
            index += 1

        label_list = source_names + target_names
        # print(label_list)

    id_list = list(range(0, len(label_list), 1))

    # Pay attention if clusters_colors or 'orig.ident_colors' missing
    if source + '_colors' not in adata.uns.keys():
        adata.uns[source + '_colors'] = list(_get_palette(adata, source).values())
    if target + '_colors' not in adata.uns.keys():
        adata.uns[target + '_colors'] = list(_get_palette(adata, target).values())

    if type(adata.uns[source + '_colors']) == np.ndarray:
        adata.uns[source + '_colors'] = adata.uns[source + '_colors'].tolist()
    if type(adata.uns[target + '_colors']) == np.ndarray:
        adata.uns[target + '_colors'] = adata.uns[target + '_colors'].tolist()

    colors = adata.uns[source + '_colors'] + adata.uns[target + '_colors']
    for number in id_list:
        tmp_list = [number, label_list[number], colors[number]]
        nodes.append(tmp_list)

    ###### LINKS ######
    key_infos_values = key_infos.values.tolist()
    key_infos_index = key_infos.index.tolist()

    # make the link df
    links = [['Source', 'Target', 'Value', 'Link Color']]
    index_target = len(label_list) - len(key_infos.index.tolist())
    for index_value in key_infos_values:
        index_source = 0
        for count in index_value:
            tmp_list = [index_source, index_target, count, colors[index_source]]
            index_source += 1
            links.append(tmp_list)
        index_target += 1

    # Retrieve headers and build dataframes
    nodes_headers = nodes.pop(0)
    links_headers = links.pop(0)
    df_nodes = pd.DataFrame(nodes, columns=nodes_headers)
    df_links = pd.DataFrame(links, columns=links_headers)

    return df_nodes, df_links


def __plot_sankey(df_nodes,
                  df_links,
                  title="Draw Sankey Diagram from dataframes",
                  save=None, scale=1):
    """
    """
    # Sankey plot setup
    data_trace = dict(
        type='sankey',
        domain=dict(
            x=[0, 1],
            y=[0, 1]
        ),
        orientation="h",
        valueformat=".0f",
        node=dict(
            pad=10,
            # thickness = 30,
            line=dict(
                color="black",
                width=0
            ),
            label=df_nodes['Label'].dropna(axis=0, how='any'),
            color=df_nodes['Color']
        ),
        link=dict(
            source=df_links['Source'].dropna(axis=0, how='any'),
            target=df_links['Target'].dropna(axis=0, how='any'),
            value=df_links['Value'].dropna(axis=0, how='any'),
            color=df_links['Link Color'].dropna(axis=0, how='any'),
        )
    )

    layout = dict(
        title=title,
        height=772,
        font=dict(
            size=10), )

    fig = dict(data=[data_trace], layout=layout)
    # fig.savefig('test.png')
    iplot(fig, validate=False)
    if save:
        pio.write_image(fig, save, width=700, height=775, scale=scale)
