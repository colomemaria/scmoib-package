import anndata as ad
import episcanpy as epi
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import plotly.io as pio
from scanpy.plotting._tools.scatterplots import _get_palette


def river_plot_2_omics(adata,
                       source,
                       target,
                       omic,
                       cell_number=False,
                       title='River plot (Sankey Diagram)',
                       save=None,
                       scale=4):
    """
    cell number count doesn't work in this function yet. 
    """

    omics = list(set(adata.obs[omic]))
    # print(omics[1], ' cells on the left')
    adata_rna = adata[adata.obs[omic] == omics[1], :].copy()
    # print(omics[0], ' cells on the right')
    adata_atac = adata[adata.obs[omic] == omics[0], :].copy()

    df_nodes_rna, df_links_rna = __tool_sankey(adata_rna,
                                               source=source,
                                               target=target,
                                               cell_number=False)

    key_infos = pd.crosstab(adata_atac.obs[target], adata_atac.obs[source])
    # key_infos

    label_list = key_infos.columns.tolist()

    id_list = list(range(df_nodes_rna.shape[0], df_nodes_rna.shape[0] + len(label_list), 1))
    nodes = [['ID', 'Label', 'Color']]

    if source + '_colors' not in adata_atac.uns.keys():
        adata_atac.uns[source + '_colors'] = list(_get_palette(adata_atac, source).values())
    if type(adata_atac.uns[source + '_colors']) == np.ndarray:
        adata_atac.uns[source + '_colors'] = adata_atac.uns[source + '_colors'].tolist()
    colors = adata.uns[source + '_colors'] + adata.uns[target + '_colors']
    index = 0
    for number in id_list:
        tmp_list = [number, label_list[index], colors[index]]
        nodes.append(tmp_list)
        index += 1

    ### merge atac nodes to rna nodes
    nodes_headers = nodes.pop(0)
    df_nodes = pd.DataFrame(nodes, columns=nodes_headers)
    df_nodes = df_nodes_rna.append(df_nodes)
    index_target = df_nodes_rna.shape[0] - len(nodes)
    # print(index_target)
    del nodes

    ### add cell number
    if cell_number == True:
        df_nodes.index = df_nodes['ID']

        key_infos = pd.crosstab(adata_atac.obs[target], adata_atac.obs[source], margins=True)
        atac_cell_numbers_source = ['(n=' + str(x) + ')' for x in key_infos.loc['All'].tolist()[:-1]]
        atac_cell_numbers_target = key_infos['All'].tolist()[:-1]
        key_infos = pd.crosstab(adata_rna.obs[target], adata_rna.obs[source], margins=True)
        rna_cell_numbers_source = ['(n=' + str(x) + ')' for x in key_infos.loc['All'].tolist()[:-1]]
        rna_cell_numbers_target = key_infos['All'].tolist()[:-1]

        rna_cell_numbers_target = [": ".join([str(omics[1]), str(x)]) for x in rna_cell_numbers_target]
        atac_cell_numbers_target = [": ".join([str(omics[0]), str(x)]) for x in atac_cell_numbers_target]
        target_cell_numbers = []
        index = 0
        for rna_count in rna_cell_numbers_target:
            target_cell_numbers.append('(' + ' & '.join([str(rna_count), str(atac_cell_numbers_target[index])]) + ')')
            index += 1

        total_count = rna_cell_numbers_source + target_cell_numbers + atac_cell_numbers_source

        new_label = []
        index = 0
        for n_cells in total_count:
            new_label.append(' '.join([str(df_nodes['Label'][index]), str(n_cells)]))
            index += 1
        df_nodes['Label'] = new_label

    ###### LINKS ######
    key_infos_values = key_infos.values.tolist()
    key_infos_index = key_infos.index.tolist()

    # make the link df
    links = [['Source', 'Target', 'Value', 'Link Color']]
    # index_target = len(label_list)-len(key_infos.index.tolist())
    # print(key_infos)
    for index_value in key_infos_values:
        index_source = df_nodes_rna.shape[0]
        index_color = 0
        for count in index_value:
            tmp_list = [index_source, index_target, count, colors[index_color]]
            index_source += 1
            index_color += 1
            links.append(tmp_list)
        index_target += 1

    ### merge atac links to rna links
    links_headers = links.pop(0)
    df_links = pd.DataFrame(links, columns=links_headers)
    tmp_var = df_links['Source'].tolist()
    tmp_var2 = df_links['Target'].tolist()
    df_links['Source'] = tmp_var2
    df_links['Target'] = tmp_var
    del tmp_var, tmp_var2

    df_links = df_links_rna.append(df_links)

    new_title = title + '\n (' + omics[1] + ' cells on the left & ' + omics[0] + ' cells on the right)'
    __plot_sankey(df_nodes, df_links,
                  title=new_title,
                  save=save,
                  scale=4)


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
