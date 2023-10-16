import numpy as np
import pandas as pd


def getDistanceFromPoint(
    adata,
    point,
    x="X_centroid",
    y="Y_centroid",
    region_col="region",
    region_subset=None,
    metric="angular",
    name="distance",
    inplace=True,
    binned=False,
    binsize=10,
):
    """
    Get the distance of each cell from a point.
    :param adata: Anndata object
    :param point: iterable coordinate of a point in (x, y) to calculate distance from
    :param x: Name of the column containing the x coordinate. Default is "X_centroid".
    :param y: Name of the column containing the y coordinate. Default is "Y_centroid".
    :param region_col: Name of the column containing the region. Default is "region".
    :param region_subset: List of regions to consider. If None, consider all cells.
    :param metric: metric to use for distance calculation.
        Metric can be "angular" or "euclidean". Default is "angular".
    :param name: Name of the column to store the distance in. Default is "distance".
    :param inplace: If True, add the distance column to adata.obs. If False, return a copy
    :param binned: If True, bin the distances into bins of size binsize.
    :param binsize: Size of the bins to use for binning. Default is 10.
    :return: If inplace is False, return a copy of adata with the distance column added
    """
    if not inplace:
        adata = adata.copy()
    adata.obs[name] = np.nan
    if region_subset is None:
        region_subset = adata.obs[region_col].unique()
    region_obs = adata.obs.loc[adata.obs[region_col].isin(region_subset), :]
    diff_x = region_obs[x] - point[0]
    diff_y = region_obs[y] - point[1]
    if metric == "euclidean":
        adata.obs.loc[adata.obs[region_col].isin(region_subset), name] = np.sqrt(
            diff_x**2 + diff_y**2
        )
    elif metric == "angular":
        adata.obs.loc[adata.obs[region_col].isin(region_subset), name] = (
            np.arctan2(diff_y, diff_x) * 180 / np.pi + 180
        )
    else:
        raise ValueError("Invalid metric. Metric must be 'angular' or 'euclidean'.")
    if binned:
        bins = adata.obs[name] - adata.obs[name] % binsize
        bins = pd.to_numeric(bins, errors="coerce").astype("Int64")
        adata.obs[name + "_binned"] = (
            "[" + bins.astype(str) + ", " + (bins + binsize).astype(str) + ")"
        )
        adata.obs.loc[
            adata.obs[name + "_binned"].isin(["[nan, nan)", "[<NA>, <NA>)"]),
            name + "_binned",
        ] = np.nan
        adata.obs.sort_values(by=name, inplace=True)
        categories = adata.obs[name + "_binned"].dropna().unique()
        adata.obs[name + "_binned"] = pd.Categorical(
            adata.obs[name + "_binned"], ordered=True, categories=categories
        )
    if not inplace:
        return adata
