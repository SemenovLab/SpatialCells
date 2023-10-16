def setGate(adata, target_marker, gate, marker_suffix="_positive", debug=True):
    """
    A new column will be added for the marker.
    E.g., for a marker named by SOX10, SOX10_positive will be added, which has two values: True or False.
    :param adata: AnnData object
    :param target_marker: the marker to be gated
    :param gate: the gate value
    :param marker_suffix: the suffix of the new column name
    :param debug: if True, print the value counts of the new column
    :return: None
    """

    newkey = target_marker + marker_suffix
    adata.obs[newkey] = (
        (adata[list(adata.obs_names), [target_marker]].X >= gate).flatten().tolist()
    )
    if debug:
        print(adata.obs[newkey].value_counts())
