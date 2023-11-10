def getRegionComposition(adata, phenotype_col, regions=None, regioncol="region"):
    """Get the cell type composition of a region.

    :param adata: Anndata object
    :param phenotype_col: list of columns containing the cell type markers
    :param regions: List of regions to consider. If None, consider all cells.
    :param regioncol: Column containing the region information
    :returns: A dataframe containing the cell type composition of the region
    """
    if regions is None:
        regions = adata.obs[regioncol].unique().tolist()
    elif not isinstance(regions, list):
        regions = [regions]
    if not isinstance(phenotype_col, list):
        phenotype_col = [phenotype_col]
    region_obs = adata.obs.loc[adata.obs[regioncol].isin(regions), :]
    region_composition = region_obs.value_counts(phenotype_col)
    region_composition = region_composition.reset_index().rename(
        columns={region_composition.name: "cell_count"}
    )
    if len(phenotype_col) == 1:
        region_composition["composition"] = (
            region_composition["cell_count"] / region_obs.shape[0]
        )
    else:
        counts = (
            region_composition.groupby(phenotype_col[:-1])["cell_count"]
            .sum()
            .reset_index()
            .rename(columns={"cell_count": "total_count"})
        )
        region_composition = region_composition.merge(
            counts, how="left", on=phenotype_col[:-1]
        )
        region_composition["composition"] = (
            region_composition["cell_count"] / region_composition["total_count"]
        )
        region_composition.drop(columns=["total_count"], inplace=True)
    return region_composition
