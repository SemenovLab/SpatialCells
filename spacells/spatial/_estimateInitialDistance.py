'''
@authors: Guihong Wan
@date: July 22, 2023
'''

import numpy as np
import scanpy as sc
from sklearn.cluster import AgglomerativeClustering

def estimateInitialDistance(adata, markders_of_interest, sampling_ratio=1.0):
    '''
    Use hierarchical clustering to get the checkpoints to estimate the distance parameter 
    for density-based clustering algorithms, e.g., DBSCAN.
    '''
    
    # get the subset data
    assert len(markders_of_interest) > 0, "markders_of_interest is empty?"
    condition = False
    for target_marker in markders_of_interest:
        try:
            newkey = target_marker+"_b"
            condition = condition | (adata.obs[newkey])
        except:
            print("please setGate for", target_marker)
            return None
    adata_tmp = adata[condition]
    
    assert (sampling_ratio >= 0) and (sampling_ratio <= 1), "sampling_ratio should be [0, 1]"
    if sampling_ratio < 1:
        sc.pp.subsample(adata_tmp, fraction=sampling_ratio, random_state=42)
    X = adata_tmp.obs[["X_centroid", "Y_centroid"]].to_numpy()
    
    # compute the distances
    model = AgglomerativeClustering(distance_threshold=0, 
        n_clusters=None,compute_distances=True)
    model = model.fit(X)
    print("Computing distances...")
    nsamples = X.shape[0]
    newnodes = []
    layer = 0
    distances = []
    for i, item in enumerate(model.children_):
        # print("i="+str(i)+":", item[0], item[1],"d:", model.distances_[i], "-->",  nsamples+i)
        newnodes.append(nsamples+i)
        if (item[0] in newnodes) and (item[1] in newnodes) :
            # print("-------",layer, model.distances_[i])
            newnodes = []
            layer += 1
            distances.append(model.distances_[i])
    # print(min(model.distances_), max(model.distances_), np.median(model.distances_), np.mean(model.distances_))
    return distances


