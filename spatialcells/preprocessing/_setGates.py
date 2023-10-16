import pandas as pd


def setGates(adata, gate_file, marker_suffix="_positive", debug=True):
    """
    A new column will be added for a marker.
    E.g., for SOX10, SOX10_positive will be added, which has two values: True or False.
    :param adata: AnnData object
    :param gate_file: the file containing the gate values for each marker
    :param marker_suffix: the suffix of the new column name
    :param debug: if True, print the value counts of the new column
    """

    manual_gate = pd.read_csv(gate_file)

    # Some markers may not be in adata
    markers_gated = []
    for m in list(adata.var_names):
        if m in list(manual_gate["markers"]):
            markers_gated.append(m)

    # Set gates
    for target_marker in markers_gated:
        val = manual_gate.loc[manual_gate["markers"] == target_marker, "gate"].values[0]
        newkey = target_marker + "_positive"
        adata.obs[newkey] = (
            (adata[list(adata.obs_names), [target_marker]].X >= val).flatten().tolist()
        )
        if debug:
            print("gate:", val, adata.obs[newkey].value_counts())
