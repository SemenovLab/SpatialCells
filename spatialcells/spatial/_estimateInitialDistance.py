"""
@authors: Guihong Wan
@date: July 22, 2023
"""

import scanpy as sc
from sklearn.cluster import AgglomerativeClustering


def estimateInitialDistance(adata, markers_of_interest, sampling_ratio=0.1):
    """
    Use hierarchical clustering to get the checkpoints to estimate the distance parameter
    for density-based clustering algorithms, e.g., DBSCAN.

    :param adata: the anndata object
    :param markers_of_interest: the list of marker names to subset the data
    :param sampling_ratio: the sampling ratio to subsample the data
    :returns: the list of distances

    """

    # get the subset data
    assert len(markers_of_interest) > 0, "markers_of_interest is empty?"
    condition = False
    for target_marker in markers_of_interest:
        try:
            condition = condition | (adata.obs[target_marker])
        except:
            print(target_marker, "is not found.")
            return None
    adata_tmp = adata[condition]

    assert (sampling_ratio >= 0) and (
        sampling_ratio <= 1
    ), "sampling_ratio should be [0, 1]"
    if sampling_ratio < 1:
        sc.pp.subsample(adata_tmp, fraction=sampling_ratio, random_state=42)
    X = adata_tmp.obs[["X_centroid", "Y_centroid"]].to_numpy()

    # compute the distances
    model = AgglomerativeClustering(
        distance_threshold=0, n_clusters=None, compute_distances=True
    )
    model = model.fit(X)
    print("Computing distances...")
    nsamples = X.shape[0]
    newnodes = []
    layer = 0
    distances = []
    for i, item in enumerate(model.children_):
        # print("i="+str(i)+":", item[0], item[1],"d:", model.distances_[i], "-->",  nsamples+i)
        newnodes.append(nsamples + i)
        if (item[0] in newnodes) and (item[1] in newnodes):
            # print("-------",layer, model.distances_[i])
            newnodes = []
            layer += 1
            distances.append(model.distances_[i])
    # print(min(model.distances_), max(model.distances_), np.median(model.distances_), np.mean(model.distances_))
    return distances
