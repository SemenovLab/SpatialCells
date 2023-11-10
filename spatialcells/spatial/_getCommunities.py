import numpy as np
from sklearn.cluster import DBSCAN


def getCommunities(
    adata,
    markers_of_interest,
    eps,
    min_samples=20,
    newcolumn="COI_community",
    core_only=False,
):
    """
    Get the communities of interest (COI) using DBSCAN

    :param adata: the anndata object
    :param markers_of_interest: the list of marker names to subset the data
    :param eps: the eps parameter for DBSCAN
    :param min_samples: Minimum number of samples in each community
    :param newcolumn: the column name of the community
    :returns: the communities of interest (COI)
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
    if core_only:
        core_samples_mask = np.zeros_like(labels, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels[~core_samples_mask] = -3
    adata.obs.loc[condition, newcolumn] = labels

    labels_sorted = []  # a list of (number of cells, cluster index)
    new_n_clusters_ = 0
    # only count clusters of interest, i.e.,
    # excluding noise (-1), not assigned (-2), and non-core (-3)
    n_clusters_ = len(set(labels) - {-1, -2, -3})
    for i in range(n_clusters_):
        idx = labels == i
        npoints = sum(idx)
        if npoints > 0:
            labels_sorted.append((npoints, i))
            new_n_clusters_ = new_n_clusters_ + 1
    labels_sorted = sorted(labels_sorted, reverse=True)

    return labels_sorted, db
