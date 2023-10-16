import numpy as np


def getMPI(
    adata,
    prolif_markers,
    arrest_markers,
    thresh_prolif=0.5,
    thresh_arrest=0.5,
    use_obs=False,
    use_layer=None,
    col_name="MPI",
    inplace=True,
):
    """
    Get MPI from a list of markers and thresholds, adapted from Gaglia et al. 2022
    https://doi.org/10.1038/s41556-022-00860-9
    MPI = :
     -1 if max(arrest_markers) > thresh_arrest
      1 else if max(prolif_markers) > thresh_prolif
      0 otherwise

    :param adata: AnnData object
    :param prolif_marker: List of proliferation markers
    :param arrest_markers: List of arrest markers
    :param thresh_prolif: Threshold for proliferation. Default is 0.5
    :param thresh_arrest: Threshold for arrest, which should be set
        based on the expression levels of KI67 marker. Default is 0.5
    :param use_obs: If True, use adata.obs[use_obs] to get the markers.
        Overrides use_layer. If use_obs==False and use_layer is None, use adata.X
    :param use_layer: Layer to use for the analysis.
        If use_obs==False and use_layer is None, use adata.X
    :param col_name: Name of the column to add to adata.obs
    :param inplace: If True, add the column to adata.obs.
        If False, return a copy of adata with the column added
    :return: None, adds a column to adata.obs
    """
    if use_obs:
        adata_df = adata.obs
        marker_names = adata.obs.columns
    elif use_layer:
        adata_df = adata.to_df(layer=use_layer)
        marker_names = adata.var_names
    else:
        adata_df = adata.to_df()
        marker_names = adata.var_names
    prolif_markers = [x for x in prolif_markers if x in marker_names]
    arrest_markers = [x for x in arrest_markers if x in marker_names]
    prolif = adata_df[prolif_markers].to_numpy()
    arrest = adata_df[arrest_markers].to_numpy()

    MPI = np.zeros(adata.n_obs)
    if prolif.shape[1] > 0:
        prolif_max = np.max(prolif, axis=1)
        MPI[prolif_max > thresh_prolif] = 1
    else:
        print("No proliferation markers found. Skipping.")

    if arrest.shape[1] > 0:
        arrest_max = np.max(arrest, axis=1)
        MPI[arrest_max > thresh_arrest] = -1
    else:
        print("No arrest markers found. Skipping.")
    if inplace:
        adata.obs[col_name] = MPI
    else:
        bdata = adata.copy()
        bdata.obs[col_name] = MPI
        return bdata
