import numpy as np


def setGate(adata, target_marker, gate, debug=True):
    '''

    @adata: AnnData
    @marker:
    @gate: the gate value for the given marker
    '''
    newkey = target_marker+"_b"
    adata.obs[newkey] = (adata[list(adata.obs_names), [target_marker]].X >= gate).flatten().tolist()
    if debug: print(adata.obs[newkey].value_counts())

