import numpy as np
import math
from ._utils import *

def getBoundary(anndata, communityColumn, communityIndexList, alpha=400, nedges_min=50, nedges_out_min=None):
    """
    Get the boundary of the community
    @anndata: the anndata object
    @communityIndexList: the list of community indexes
    @alpha: the alpha parameter for alpha shape
    @nedges_min: the minimum number of edges for a community
    @nedges_out_min: the minimum number of edges for a community to be considered as "out"

    Return: the boundary of the community
    """
    xy = anndata.obs[anndata.obs[communityColumn].isin(communityIndexList)][
        ["X_centroid", "Y_centroid"]
    ].to_numpy()
    print("xy shape:", xy.shape)
    edge_components = getAlphaShapes(xy, alpha)
    print("edge_components:", len(edge_components), [len(i) for i in edge_components], edge_components[0].shape)
    grouped_components = groupRemoveEdgeComponents(edge_components, nedges_min, nedges_out_min)
    return grouped_components

