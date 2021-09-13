from scanpy.plotting._tools.scatterplots import _get_palette
from scanpy._utils import sanitize_anndata, _doc_params, Empty, _empty
from scanpy.plotting._tools.scatterplots import _get_data_points, _get_color_source_vector, _utils
from scanpy.plotting._tools.scatterplots import _color_vector, settings, _get_vmin_vmax, _basis2name
import scanpy.plotting._utils
import scanpy._utils
import logging
from scanpy._utils import NeighborsView
from scanpy.plotting._utils import _get_basis
from scanpy.plotting._tools.scatterplots import _add_categorical_legend
import warnings
import itertools
from packaging.version import parse
import networkx as nx

import collections.abc as cabc
from copy import copy
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple, Callable

import numpy as np
import pandas as pd
from anndata import AnnData
from cycler import Cycler
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from pandas.api.types import is_categorical_dtype
import matplotlib
from matplotlib import pyplot as pl, colors
from matplotlib.cm import get_cmap
from matplotlib import rcParams
from matplotlib import patheffects
from matplotlib.colors import Colormap
from functools import partial

try:
    from typing import Literal
except ImportError:
    try:
        from typing_extensions import Literal
    except ImportError:

        class LiteralMeta(type):
            def __getitem__(cls, values):
                if not isinstance(values, tuple):
                    values = (values,)
                return type('Literal_', (Literal,), dict(__args__=values))


        class Literal(metaclass=LiteralMeta):
            pass

ColorLike = Union[str, Tuple[float, ...]]
_IGraphLayout = Literal['fa', 'fr', 'rt', 'rt_circular', 'drl', 'eq_tree', ...]
_FontWeight = Literal['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']
_FontSize = Literal[
    'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'
]

VMinMax = Union[str, float, Callable[[Sequence[float]], float]]


def umap_barcodes(adata, bc_list1, bc_list2, basis='umap', color=None, edges=True, edges_width=0.1):
    bc_list = list(adata.obs.index)
    edge_list = list(zip([bc_list.index(i) for i in bc_list1], [bc_list.index(j) for j in bc_list2]))
    __embedding(adata, basis=basis, color=color, edges=edges, edges_width=edges_width, edge_list=edge_list)


def __check_projection(projection):
    """Validation for projection argument."""
    if projection not in {"2d", "3d"}:
        raise ValueError(f"Projection must be '2d' or '3d', was '{projection}'.")
    if projection == "3d":
        mpl_version = parse(matplotlib.__version__)
        if mpl_version < parse("3.3.3"):
            raise ImportError(
                f"3d plotting requires matplotlib > 3.3.3. Found {matplotlib.__version__}"
            )


def __plot_edges(axs, adata, basis, edges_width, edges_color, edge_list=None, neighbors_key=None):
    if not isinstance(axs, cabc.Sequence):
        axs = [axs]

    if neighbors_key is None:
        neighbors_key = 'neighbors'
    if neighbors_key not in adata.uns:
        raise ValueError('`edges=True` requires `pp.neighbors` to be run before.')
    neighbors = NeighborsView(adata, neighbors_key)
    g = nx.Graph(neighbors['connectivities'])
    basis_key = _get_basis(adata, basis)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for ax in axs:
            edge_collection = nx.draw_networkx_edges(
                g,
                adata.obsm[basis_key],
                edgelist=edge_list,
                ax=ax,
                width=edges_width,
                edge_color=edges_color,
            )
            edge_collection.set_zorder(-2)
            edge_collection.set_rasterized(settings._vector_friendly)


def __embedding(
        adata: AnnData,
        basis: str,
        *,
        color: Union[str, Sequence[str], None] = None,
        gene_symbols: Optional[str] = None,
        use_raw: Optional[bool] = None,
        sort_order: bool = True,
        edges: bool = False,
        edges_width: float = 0.1,
        edges_color: Union[str, Sequence[float], Sequence[str]] = 'grey',
        edge_list=None,
        neighbors_key: Optional[str] = None,
        arrows: bool = False,
        arrows_kwds: Optional[Mapping[str, Any]] = None,
        groups: Optional[str] = None,
        components: Union[str, Sequence[str]] = None,
        layer: Optional[str] = None,
        projection: Literal['2d', '3d'] = '2d',
        scale_factor: Optional[float] = None,
        color_map: Union[Colormap, str, None] = None,
        cmap: Union[Colormap, str, None] = None,
        palette: Union[str, Sequence[str], Cycler, None] = None,
        na_color: ColorLike = "lightgray",
        na_in_legend: bool = True,
        size: Union[float, Sequence[float], None] = None,
        frameon: Optional[bool] = None,
        legend_fontsize: Union[int, float, _FontSize, None] = None,
        legend_fontweight: Union[int, _FontWeight] = 'bold',
        legend_loc: str = 'right margin',
        legend_fontoutline: Optional[int] = None,
        vmax: Union[VMinMax, Sequence[VMinMax], None] = None,
        vmin: Union[VMinMax, Sequence[VMinMax], None] = None,
        add_outline: Optional[bool] = False,
        outline_width: Tuple[float, float] = (0.3, 0.05),
        outline_color: Tuple[str, str] = ('black', 'white'),
        ncols: int = 4,
        hspace: float = 0.25,
        wspace: Optional[float] = None,
        title: Union[str, Sequence[str], None] = None,
        show: Optional[bool] = None,
        save: Union[bool, str, None] = None,
        ax: Optional[Axes] = None,
        return_fig: Optional[bool] = None,
        **kwargs,
) -> Union[Figure, Axes, None]:
    """\
    Scatter plot for user specified embedding basis (e.g. umap, pca, etc)

    Parameters
    ----------
    basis
        Name of the `obsm` basis to use.
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    __check_projection(projection)
    sanitize_anndata(adata)

    # Setting up color map for continuous values
    if color_map is not None:
        if cmap is not None:
            raise ValueError("Cannot specify both `color_map` and `cmap`.")
        else:
            cmap = color_map
    cmap = copy(get_cmap(cmap))
    cmap.set_bad(na_color)
    kwargs["cmap"] = cmap

    # Prevents warnings during legend creation
    na_color = colors.to_hex(na_color, keep_alpha=True)

    if size is not None:
        kwargs['s'] = size
    if 'edgecolor' not in kwargs:
        # by default turn off edge color. Otherwise, for
        # very small sizes the edge will not reduce its size
        # (https://github.com/theislab/scanpy/issues/293)
        kwargs['edgecolor'] = 'none'

    if groups:
        if isinstance(groups, str):
            groups = [groups]

    args_3d = dict(projection='3d') if projection == '3d' else {}

    # Deal with Raw
    if use_raw is None:
        # check if adata.raw is set
        use_raw = layer is None and adata.raw is not None
    if use_raw and layer is not None:
        raise ValueError(
            "Cannot use both a layer and the raw representation. Was passed:"
            f"use_raw={use_raw}, layer={layer}."
        )

    if wspace is None:
        #  try to set a wspace that is not too large or too small given the
        #  current figure size
        wspace = 0.75 / rcParams['figure.figsize'][0] + 0.02
    if adata.raw is None and use_raw:
        raise ValueError(
            "`use_raw` is set to True but AnnData object does not have raw. "
            "Please check."
        )
    # turn color into a python list
    color = [color] if isinstance(color, str) or color is None else list(color)
    if title is not None:
        # turn title into a python list if not None
        title = [title] if isinstance(title, str) else list(title)

    # get the points position and the components list
    # (only if components is not None)
    data_points, components_list = _get_data_points(
        adata, basis, projection, components, scale_factor
    )

    # Setup layout.
    # Most of the code is for the case when multiple plots are required
    # 'color' is a list of names that want to be plotted.
    # Eg. ['Gene1', 'louvain', 'Gene2'].
    # component_list is a list of components [[0,1], [1,2]]
    if (
            not isinstance(color, str)
            and isinstance(color, cabc.Sequence)
            and len(color) > 1
    ) or len(components_list) > 1:
        if ax is not None:
            raise ValueError(
                "Cannot specify `ax` when plotting multiple panels "
                "(each for a given value of 'color')."
            )
        if len(components_list) == 0:
            components_list = [None]

        # each plot needs to be its own panel
        num_panels = len(color) * len(components_list)
        fig, grid = _panel_grid(hspace, wspace, ncols, num_panels)
    else:
        if len(components_list) == 0:
            components_list = [None]
        grid = None
        if ax is None:
            fig = pl.figure()
            ax = fig.add_subplot(111, **args_3d)

    # turn vmax and vmin into a sequence
    if isinstance(vmax, str) or not isinstance(vmax, cabc.Sequence):
        vmax = [vmax]
    if isinstance(vmin, str) or not isinstance(vmin, cabc.Sequence):
        vmin = [vmin]

    if 's' in kwargs:
        size = kwargs.pop('s')

    if size is not None:
        # check if size is any type of sequence, and if so
        # set as ndarray
        import pandas.core.series

        if (
                size is not None
                and isinstance(size, (cabc.Sequence, pandas.core.series.Series, np.ndarray))
                and len(size) == adata.shape[0]
        ):
            size = np.array(size, dtype=float)
    else:
        size = 120000 / adata.shape[0]

    ###
    # make the plots
    axs = []

    idx_components = range(len(components_list))

    # use itertools.product to make a plot for each color and for each component
    # For example if color=[gene1, gene2] and components=['1,2, '2,3'].
    # The plots are: [
    #     color=gene1, components=[1,2], color=gene1, components=[2,3],
    #     color=gene2, components = [1, 2], color=gene2, components=[2,3],
    # ]
    for count, (value_to_plot, component_idx) in enumerate(
            itertools.product(color, idx_components)
    ):
        color_source_vector = _get_color_source_vector(
            adata,
            value_to_plot,
            layer=layer,
            use_raw=use_raw,
            gene_symbols=gene_symbols,
            groups=groups,
        )
        color_vector, categorical = _color_vector(
            adata,
            value_to_plot,
            color_source_vector,
            palette=palette,
            na_color=na_color,
        )

        ### Order points
        order = slice(None)
        if sort_order is True and value_to_plot is not None and categorical is False:
            # Higher values plotted on top, null values on bottom
            order = np.argsort(-color_vector, kind="stable")[::-1]
        elif sort_order and categorical:
            # Null points go on bottom
            order = np.argsort(~pd.isnull(color_source_vector), kind="stable")
        # Set orders
        if isinstance(size, np.ndarray):
            size = np.array(size)[order]
        color_source_vector = color_source_vector[order]
        color_vector = color_vector[order]
        _data_points = data_points[component_idx][order, :]

        # if plotting multiple panels, get the ax from the grid spec
        # else use the ax value (either user given or created previously)
        if grid:
            ax = pl.subplot(grid[count], **args_3d)
            axs.append(ax)
        if not (settings._frameon if frameon is None else frameon):
            ax.axis('off')
        if title is None:
            if value_to_plot is not None:
                ax.set_title(value_to_plot)
            else:
                ax.set_title('')
        else:
            try:
                ax.set_title(title[count])
            except IndexError:
                logging.warning(
                    "The title list is shorter than the number of panels. "
                    "Using 'color' value instead for some plots."
                )
                ax.set_title(value_to_plot)

        # check vmin and vmax options
        if categorical:
            kwargs['vmin'] = kwargs['vmax'] = None
        else:
            kwargs['vmin'], kwargs['vmax'] = _get_vmin_vmax(
                vmin, vmax, count, color_vector
            )

        # make the scatter plot
        if projection == '3d':
            cax = ax.scatter(
                _data_points[:, 0],
                _data_points[:, 1],
                _data_points[:, 2],
                marker=".",
                c=color_vector,
                rasterized=settings._vector_friendly,
                **kwargs,
            )
        else:

            scatter = (
                partial(ax.scatter, s=size, plotnonfinite=True)
                if scale_factor is None
                else partial(circles, s=size, ax=ax)  # size in circles is radius
            )

            if add_outline:
                # the default outline is a black edge followed by a
                # thin white edged added around connected clusters.
                # To add an outline
                # three overlapping scatter plots are drawn:
                # First black dots with slightly larger size,
                # then, white dots a bit smaller, but still larger
                # than the final dots. Then the final dots are drawn
                # with some transparency.

                bg_width, gap_width = outline_width
                point = np.sqrt(size)
                gap_size = (point + (point * gap_width) * 2) ** 2
                bg_size = (np.sqrt(gap_size) + (point * bg_width) * 2) ** 2
                # the default black and white colors can be changes using
                # the contour_config parameter
                bg_color, gap_color = outline_color

                # remove edge from kwargs if present
                # because edge needs to be set to None
                kwargs['edgecolor'] = 'none'

                # remove alpha for outline
                alpha = kwargs.pop('alpha') if 'alpha' in kwargs else None

                ax.scatter(
                    _data_points[:, 0],
                    _data_points[:, 1],
                    s=bg_size,
                    marker=".",
                    c=bg_color,
                    rasterized=settings._vector_friendly,
                    **kwargs,
                )
                ax.scatter(
                    _data_points[:, 0],
                    _data_points[:, 1],
                    s=gap_size,
                    marker=".",
                    c=gap_color,
                    rasterized=settings._vector_friendly,
                    **kwargs,
                )
                # if user did not set alpha, set alpha to 0.7
                kwargs['alpha'] = 0.7 if alpha is None else alpha

            cax = scatter(
                _data_points[:, 0],
                _data_points[:, 1],
                marker=".",
                c=color_vector,
                rasterized=settings._vector_friendly,
                **kwargs,
            )

        # remove y and x ticks
        ax.set_yticks([])
        ax.set_xticks([])
        if projection == '3d':
            ax.set_zticks([])

        # set default axis_labels
        name = _basis2name(basis)
        if components is not None:
            axis_labels = [name + str(x + 1) for x in components_list[component_idx]]
        elif projection == '3d':
            axis_labels = [name + str(x + 1) for x in range(3)]
        else:
            axis_labels = [name + str(x + 1) for x in range(2)]

        ax.set_xlabel(axis_labels[0])
        ax.set_ylabel(axis_labels[1])
        if projection == '3d':
            # shift the label closer to the axis
            ax.set_zlabel(axis_labels[2], labelpad=-7)
        ax.autoscale_view()

        if edges:
            __plot_edges(ax, adata, basis, edges_width, edges_color, edge_list, neighbors_key)
        if arrows:
            _utils.plot_arrows(ax, adata, basis, arrows_kwds)

        if value_to_plot is None:
            # if only dots were plotted without an associated value
            # there is not need to plot a legend or a colorbar
            continue

        if legend_fontoutline is not None:
            path_effect = [
                patheffects.withStroke(linewidth=legend_fontoutline, foreground='w')
            ]
        else:
            path_effect = None

        # Adding legends
        if categorical:
            _add_categorical_legend(
                ax,
                color_source_vector,
                palette=_get_palette(adata, value_to_plot),
                scatter_array=_data_points,
                legend_loc=legend_loc,
                legend_fontweight=legend_fontweight,
                legend_fontsize=legend_fontsize,
                legend_fontoutline=path_effect,
                na_color=na_color,
                na_in_legend=na_in_legend,
                multi_panel=bool(grid),
            )
        else:
            # TODO: na_in_legend should have some effect here
            pl.colorbar(cax, ax=ax, pad=0.01, fraction=0.08, aspect=30)

    if return_fig is True:
        return fig
    axs = axs if grid else ax
    _utils.savefig_or_show(basis, show=show, save=save)
    if show is False:
        return axs
