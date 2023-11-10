from ._getRegionArea import getRegionArea


def getRegionDensity(
    adata,
    boundary,
    region_col="region",
    region_subset=None,
    phenotype_col=[],
    exclude_holes=True,
):
    """Get the density of cells in a region defined by a list of Polygon objects.
    If phenotype_col is empty, return the total density.
    If exclude_holes is True, the area of the holes in the region is
    subtracted from the overall area of the region for density calculation.

    :param adata: Anndata object
    :param boundary: A multiPolygon object defining the boundary of the region
    :param region_col: Name of the column containing the region. Default is "region".
    :param region_subset: List of regions to consider. If None, consider all cells.
    :param phenotype_col: A list of columns to stratify the density by.
        If empty, return the total density.
    :param exclude_holes: whether to exclude the holes in the region
    :returns: density of cells in the region as a pandas Series
        stratified by phenotype_col
    """
    area = getRegionArea(boundary, exclude_holes)
    if region_subset is None:
        region_subset = adata.obs[region_col].unique()
    else:
        region_subset = region_subset
    adata = adata.copy()
    adata.obs["All"] = "All cells"
    if len(phenotype_col) == 0:
        phenotype_col = ["All"]
    return (
        adata.obs[adata.obs[region_col].isin(region_subset)].value_counts(phenotype_col)
        / area
    )
