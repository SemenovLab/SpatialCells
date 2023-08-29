import numpy as np
import math
from sklearn.cluster import DBSCAN


def getCommunities(
    adata, markers_of_interest, eps, min_samples=20, newcolumn="COI_community"
):
    """
    Get the communities of interest (COI) using DBSCAN
    :param adata: the anndata object
    :param markers_of_interest: the list of marker names to subset the data
    :param eps: the eps parameter for DBSCAN
    :param min_samples: Minimum number of samples in each community
    :param newcolumn: the column name of the community
    :return: the communities of interest (COI)
    """
    # get cells of interest (COI)
    assert len(markers_of_interest) > 0, "markers_of_interest is empty?"
    condition = False
    for target_marker in markers_of_interest:
        try:
            condition = condition | (adata.obs[target_marker])
        except:
            print(target_marker, "is not found.")
            return None
    adata_tmp = adata[condition]

    # call DBSCAN
    X = adata_tmp.obs[["X_centroid", "Y_centroid"]].to_numpy()
    db = DBSCAN(eps=eps, min_samples=min_samples, algorithm="kd_tree").fit(X)

    # assign labels
    labels = db.labels_
    adata.obs[newcolumn] = -2
    adata.obs.loc[condition, newcolumn] = labels

    labels_sorted = []  # a list of (number of cells, cluster index)
    new_n_clusters_ = 0
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    for i in range(n_clusters_):
        idx = labels == i
        npoints = sum(idx)
        if npoints > 0:
            labels_sorted.append((npoints, i))
            new_n_clusters_ = new_n_clusters_ + 1
    labels_sorted = sorted(labels_sorted, reverse=True)
    return labels_sorted, db
