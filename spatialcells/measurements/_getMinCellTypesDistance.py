from sklearn.metrics import pairwise_distances_argmin_min


def getMinCellTypesDistance(adata1, adata2):
    """Return the minimum distance between cell types of two AnnData objects.

    :param adata1: Anndata object
    :param adata2: Anndata object
    :returns: minimum distance between cell types
    """
    # Get the cell types of the two AnnData objects
    cell_types1 = adata1.obs[["X_centroid", "Y_centroid"]]
    cell_types2 = adata2.obs[["X_centroid", "Y_centroid"]]

    # Get the minimum distance between cell types
    _, min_distance = pairwise_distances_argmin_min(cell_types1, cell_types2)

    return min_distance
