import anndata as ad
import numpy as np
import pandas as pd

from ._getRegionComposition import getRegionComposition


def getSlidingWindowsComposition(
    adata,
    window_size,
    step_size,
    phenotype_col,
    region_col="region",
    region_subset=None,
    min_cells=0,
):
    """Get Sliding window cell composition for cells in region subset.

    :param adata: Anndata object
    :param window_size: Size of the sliding window
    :param step_size: Size of the step
    :param phenotype_col: list of columns containing the cell type markers,
        for cell type composition
    :param region_col: Column containing the region information
    :param region_subset: List of regions to consider. If None, consider all cells.
    :param min_cells: Minimum number of cells in a window to consider it
    :returns: A dataframe containing the cell type composition of the region in each window
    """
    if region_subset is None:
        obs = adata.obs
    else:
        obs = adata.obs[adata.obs[region_col].isin(region_subset)]
    if obs.empty:
        return pd.DataFrame()

    minx = int(obs["X_centroid"].min())
    miny = int(obs["Y_centroid"].min())
    maxx = int(obs["X_centroid"].max())
    maxy = int(obs["Y_centroid"].max())

    def composition_for(window_obs, x, y):
        comp = getRegionComposition(
            ad.AnnData(obs=window_obs), phenotype_col, regioncol=region_col
        )
        comp["X_start"] = x
        comp["Y_start"] = y
        return comp

    rows = []
    if window_size == step_size:
        # speed up: each cell falls into exactly one window, so bucket once
        # and group, rather than iterating all (x, y) pairs
        x_bin = ((obs["X_centroid"] - minx) // step_size * step_size + minx).astype(int)
        y_bin = ((obs["Y_centroid"] - miny) // step_size * step_size + miny).astype(int)
        for (x, y), window_obs in obs.assign(_X_bin=x_bin, _Y_bin=y_bin).groupby(
            ["_X_bin", "_Y_bin"], observed=True
        ):
            if window_obs.shape[0] > min_cells:
                rows.append(composition_for(window_obs, int(x), int(y)))
    else:
        for x in range(minx, maxx + window_size, step_size):
            for y in range(miny, maxy + window_size, step_size):
                window_obs = obs[
                    (obs["X_centroid"] >= x)
                    & (obs["X_centroid"] < x + window_size)
                    & (obs["Y_centroid"] >= y)
                    & (obs["Y_centroid"] < y + window_size)
                ]
                if window_obs.shape[0] > min_cells:
                    rows.append(composition_for(window_obs, x, y))

    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows)
    out["window_size"] = window_size
    out["step_size"] = step_size
    return out


def get_comp_mask(df, pheno_col, pheno_vals, step_size):
    """
    Get a mask of the composition of the region in each window

    :param df: A dataframe containing the cell type composition of pheno_vals in each window
    :param pheno_col: Column containing the cell type information
    :param pheno_vals: List of cell types to consider
    :param step_size: Size of the step
    :return: A np array mask of the composition of the region in each window
    """
    maxx, maxy = df["X_start"].max() + step_size, df["Y_start"].max() + step_size
    mask = np.zeros((maxy + 2000, maxx + 2000))
    df1 = df[df[pheno_col].isin(pheno_vals)]
    for i in range(len(df1)):
        x = int(df1.iloc[i]["X_start"])
        y = int(df1.iloc[i]["Y_start"])
        mask[y : y + step_size, x : x + step_size] = df1.iloc[i]["composition"]
    return mask
