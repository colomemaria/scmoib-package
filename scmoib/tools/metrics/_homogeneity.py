import episcanpy as epi


def homogeneity(adata, group1, group2):
    return epi.tl.homogeneity(adata, group1, group2)
