from sklearn.decomposition import PCA


def compute_PCs(X):
    pca = PCA(n_components=2, svd_solver='full').fit(X)
    return pca

def compute_PC_distances(pca, X):
    """[summary]

    Args:
        points (numpy.ndarray): points

    Returns:
        (float, float): explained variances on the two PCA components.
    """
    W = pca.components_.T@X.T
    d1 = abs(min(W[0])) + abs(max(W[0]))
    d2 = abs(min(W[1])) + abs(max(W[1]))
    return (d1, d2)

    
    