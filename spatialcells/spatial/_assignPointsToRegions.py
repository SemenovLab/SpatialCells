import numpy as np
import pandas as pd
from shapely.geometry import Polygon, Point
from ._utils import *


def isInPolygon(polygons, point):
    """
    Check if a point is in any of the polygons
    :param polygons: List of Polygon objects
    :param point: iterable in the form of (x, y)
    :return: True if the point is in any of the polygons, False otherwise
    """
    for polygon in polygons:
        if polygon.contains(Point(point)):
            return True
    return False


def assignPointsToRegion(
    anndata, boundaries_component, assigncolumn="region", target="target", donelist=[]
):
    """
    Assign points to a region based on the boundaries. Points that are already assigned
    to a region in donelist will not be reassigned.
    :param anndata: Anndata object
    :param boundaries_component: Boundaries of the region
    :param assigncolumn: Column name for the region assignment
    :param target: Region name to assign to the points
    """
    assign_subset = anndata.obs[~anndata.obs[assigncolumn].isin(donelist)][
        ["X_centroid", "Y_centroid", assigncolumn]
    ]
    assign_coords = assign_subset[["X_centroid", "Y_centroid"]].to_numpy()
    assign_region = assign_subset[assigncolumn].to_list()
    polygons = getPolygons(boundaries_component)
    for idx in range(len(assign_coords)):
        if isInPolygon(polygons, assign_coords[idx, :]):
            assign_region[idx] = target
    anndata.obs.loc[
        ~anndata.obs[assigncolumn].isin(donelist), assigncolumn
    ] = assign_region


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
