import episcanpy as epi


def ami(adata, group1, group2):
    return epi.tl.AMI(adata, group1, group2)