from ._accuracy_per_cell_type import accuracy_per_cell_type
from ._cd_ratio import cd_ratio
from ._cumulative_node import cumulative_node
from ._cumulative_node_group import cumulative_node_group
from ._dists_distr import dists_distr
from ._mean_node_per_cell_type import mean_node_per_cell_type
from ._metric_heatmap import metric_heatmap
from ._node_distr import node_distr
from ._river_plot import river_plot
from ._river_plot_2_omics import river_plot_2_omics

from ._pairwise_distance import pairwise_distance
from ._silhouette import silhouette
from ._metrics_scatterplot import metrics_scatterplot

import scanpy
if scanpy.__version__<1.8.0 :
	from ._umap_barcodes import umap_barcodes
else: 
	from ._umap_barcodes_2 import umap_barcodes
	