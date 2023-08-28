from ._utils import *


def getBoundary(
    anndata,
    communitycolumn,
    communityIndexList,
    alpha=100,
    nedges_min=20,
    nedges_out_min=20,
):
    """
    Get the boundary of the communities of interest
    :param anndata: the anndata object
    :param communitycolumn: the column name of the community
    :param communityIndexList: the list of community indexes
    :param alpha: the alpha parameter for alpha shape. Smaller alpha gives
        more detailed boundary, but may appear more jagged.
    :param nedges_min: the minimum number of edges to keep a hole inside a
        component's boundary
    :param nedges_out_min: the minimum number of edges for a component to be
        considered as an "out" region
    :return: the boundaries of components.
    """

    xy = anndata.obs[anndata.obs[communitycolumn].isin(communityIndexList)][
        ["X_centroid", "Y_centroid"]
    ].to_numpy()
    edge_components = getAlphaShapes(xy, alpha)
    grouped_components = groupRemoveEdgeComponents(
        edge_components, nedges_min, nedges_out_min
    )

    boundary_edges = getEdgesOnBoundaries(grouped_components)
    return boundary_edges, grouped_components
