import numpy as np
from sklearn.neighbors import KernelDensity


def getCellKDE(
    adata,
    regions,
    phenotype_col=None,
    phenotype_subset=[],
    bandwidth=1,
    name="kde_likelihood",
):
    """
    Get per cell log likelihood based on a kernel density estimate.
    This can be normalized by the area of the region to be comparable across regions.
    Likelihoods will be stored in adata.obs[name].
    :param adata: Anndata object
    :param regions: A list of regions to compute the density in
    :param phenotype_col: A list of columns to stratify the density by.
    :param phenotype_subset: A list of cell type markers to subset the data by
    """
    subset = adata.obs[
        adata.obs["region"].isin(regions)
        & adata.obs[phenotype_col].isin(phenotype_subset)
    ][["X_centroid", "Y_centroid"]]
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(subset)
    adata.obs[name] = kde.score_samples(adata.obs[["X_centroid", "Y_centroid"]])
    return kde
