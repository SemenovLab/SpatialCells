import shapely
from shapely.geometry import Point, MultiPoint
import numpy as np
import pandas as pd
from tqdm import tqdm


def getDistanceFromObject(
    adata,
    object,
    x="X_centroid",
    y="Y_centroid",
    region_col="region",
    region_subset=None,
    name="distance",
    inplace=True,
    binned=False,
    binsize=10,
):
    """
    Get the minimum euclidean distance between each cell and a shapely object.
    :param adata: Anndata object
    :param object: Shapely object to measure distance from
    :param x: Name of the column containing the x coordinate. Default is "X_centroid".
    :param y: Name of the column containing the y coordinate. Default is "Y_centroid".
    :param region_col: Name of the column containing the region. Default is "region".
    :param region_subset: List of regions to consider. If None, consider all cells.
    :param name: Name of the column to store the distance in. Default is "distance".
    :param inplace: If True, add the distance column to adata.obs. If False, return a copy
    :param binned: If True, bin the distances into bins of size binsize.
    :param binsize: Size of the bins to use for binning. Default is 10.
    :return: If inplace is False, return a copy of adata with the distance column added
    """
    if not inplace:
        adata = adata.copy()
    if not isinstance(object, shapely.geometry.base.BaseGeometry):
        raise ValueError("object must be a shapely object")
    if region_subset is None:
        region_subset = adata.obs[region_col].unique()
    region_obs = adata.obs.loc[adata.obs[region_col].isin(region_subset), :]
    points = MultiPoint(region_obs[[x, y]].to_numpy())
    dists = []
    for i, pt in tqdm(enumerate(points.geoms)):
        dists.append(Point(pt).distance(object))
    adata.obs[name] = np.nan
    adata.obs.loc[adata.obs[region_col].isin(region_subset), name] = dists
    if binned:
        bins = adata.obs[name] - adata.obs[name] % binsize
        adata.obs[name + "_binned"] = (
            "["
            + bins.astype(str).str[:-2]
            + ", "
            + (bins + binsize).astype(str).str[:-2]
            + ")"
        )
        adata.obs.loc[
            adata.obs[name + "_binned"] == "[n, n)", name + "_binned"
        ] = np.nan
        adata.obs.sort_values(by=name, inplace=True)
        categories = adata.obs[name + "_binned"].dropna().unique()
        adata.obs[name + "_binned"] = pd.Categorical(
            adata.obs[name + "_binned"], ordered=True, categories=categories
        )
    if not inplace:
        return adata
