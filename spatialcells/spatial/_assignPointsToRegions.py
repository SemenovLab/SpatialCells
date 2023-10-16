import pandas as pd
from shapely.geometry import MultiPoint
from tqdm import tqdm
from ._utils import *


def assignPointsToRegion(
    anndata, multi_polygon, assigncolumn="region", target="target", donelist=[]
):
    """
    Assign points to a region based on the boundaries. Points that are already assigned
    to a region in donelist will not be reassigned.
    :param anndata: Anndata object
    :param multi_polygon: MultiPolygon object of the region boundary
    :param assigncolumn: Column name for the region assignment
    :param target: Region name to assign to the points in the region
    :param donelist: List of region names that are already assigned and should not be reassigned
    """
    minx, miny, maxx, maxy = multi_polygon.bounds
    condition = (
        ~anndata.obs[assigncolumn].isin(donelist)
        & (anndata.obs["X_centroid"] >= minx)
        & (anndata.obs["X_centroid"] <= maxx)
        & (anndata.obs["Y_centroid"] >= miny)
        & (anndata.obs["Y_centroid"] <= maxy)
    )
    assign_subset = anndata.obs[condition][["X_centroid", "Y_centroid", assigncolumn]]
    assign_coords = assign_subset[["X_centroid", "Y_centroid"]].to_numpy()
    assign_points = MultiPoint(assign_coords)
    assign_region = assign_subset[assigncolumn].to_list()
    for i, point in tqdm(enumerate(assign_points.geoms)):
        if multi_polygon.contains(point):
            assign_region[i] = target
    anndata.obs.loc[condition, assigncolumn] = assign_region


def assignPointsToRegions(
    anndata, boundaries_list, region_names, assigncolumn="region", default="BG"
):
    """
    Assign points to regions based on the boundaries. The region assignment is
    based on the order of the boundaries, so the innermost region should be the
    first element of boundaries_list.
    :param anndata: Anndata object
    :param boundaries_list: List of boundaries
    :param region_names: List of region names. The order and length should match boundaries_list
    :param assigncolumn: Column name for the region assignment
    :param default: Default region name for points that are not assigned to any region
    """
    assert (
        len(boundaries_list) > 0
    ), "Length of boundaries_list should be greater than 0"
    assert len(boundaries_list) == len(
        region_names
    ), "Length of boundaries_list and region_names should match"
    anndata.obs[assigncolumn] = default
    anndata.obs[assigncolumn] = pd.Categorical(
        anndata.obs[assigncolumn], categories=region_names + [default], ordered=True
    )
    for idx in range(len(boundaries_list)):
        assignPointsToRegion(
            anndata,
            boundaries_list[idx],
            assigncolumn,
            region_names[idx],
            donelist=region_names[:idx],
        )
        print("Assigned points to region: " + region_names[idx])
