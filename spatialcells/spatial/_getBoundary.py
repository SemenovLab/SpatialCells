import numpy as np
import math
from shapely.geometry import Polygon, MultiPolygon
from ._utils import *

def getBoundary(anndata, 
                communitycolumn, 
                communityIndexList, 
                alpha=100, 
                nedges_min=20, 
                nedges_out_min=20):
    """
    Get the boundary of the communities
    
    @anndata: the anndata object
    @communityIndexList: the list of community indexes
    @alpha: the alpha parameter for alpha shape
    @nedges_min: the minimum number of edges for a component
    @nedges_out_min: the minimum number of edges for a component to be considered as an "out" region

    Return: the boundaries of components.
    """
    
    xy = anndata.obs[anndata.obs[communitycolumn].isin(communityIndexList)][
        ["X_centroid", "Y_centroid"]
    ].to_numpy()
    # print("xy shape:", xy.shape)
    
    edge_components = getAlphaShapes(xy, alpha)
    # return xy[edge_components]
    # return edge_components
    # print("edge_components:", len(edge_components), [len(i) for i in edge_components], edge_components[0].shape)
    
    grouped_components = groupRemoveEdgeComponents(edge_components, nedges_min, nedges_out_min)

    

    # comp_polygons = getPolygons(grouped_components)
    boundary_edges = getEdgesOnBoundaries(grouped_components)

    # return boundary_edges, grouped_components, comp_polygons
    return boundary_edges, grouped_components



