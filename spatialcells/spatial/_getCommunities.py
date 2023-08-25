import numpy as np
import math
from sklearn.cluster import DBSCAN


def getCommunities(adata, markders_of_interest, 
                   eps, min_samples = 20, 
                   newcolumn = "COI_community"):
    '''

    @adata: AnnData
    @markders_of_interest: 
    
    '''
    
    # get cells of interest (COI)
    assert len(markders_of_interest) > 0, "markders_of_interest is empty?"
    condition = False
    for target_marker in markders_of_interest:
        try:
            condition = condition | (adata.obs[target_marker])
        except:
            print(target_marker, "is not found.")
            return None
    adata_tmp = adata[condition]
    
    # call DBSCAN
    X = adata_tmp.obs[["X_centroid", "Y_centroid"]].to_numpy()
    db = DBSCAN(eps = eps, min_samples = min_samples, algorithm="kd_tree").fit(X)
    
    # assign labels 
    labels = db.labels_
    adata.obs[newcolumn] = -2
    adata.obs.loc[condition, newcolumn] = labels
    
    labels_sorted = [] # a list of (number of cells, cluster index)
    new_n_clusters_ = 0
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    for i in range(n_clusters_):
        idx = labels==i
        npoints = sum(idx)
        if npoints > 0: 
            labels_sorted.append((npoints,i))
            new_n_clusters_ = new_n_clusters_+1
    labels_sorted = sorted(labels_sorted, reverse=True)
    return labels_sorted, db

