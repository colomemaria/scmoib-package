import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def dists_distr(adata, title=None, bins='auto', figsize=(10, 5)):
    """
    Historgram plotting the weighted distance between matching barcodes. 
    Parameters
    ----------
    adata: AnnData object
    
    title: str, optional (Defailt: None)
    
    bins: str, number, vector, or a pair of such values
        Generic bin parameter that can be the name of a reference rule, the number of bins, or the breaks of the bins. 
        Passed to numpy.histogram_bin_edges() (https://seaborn.pydata.org/generated/seaborn.histplot.html).
        
    figsize: (float, float)
    """

    fig, ax = plt.subplots(figsize=figsize)
    sns.histplot(data=adata.uns['node_metrics']['dists'], ax=ax)
    ax.set(xlabel='Weighted distance between barcodes', title='Weighted distance distribution');