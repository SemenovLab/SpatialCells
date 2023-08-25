import numpy as np


def setGate(adata, target_marker, gate, debug=True):
    '''
    A new column will be added for the marker. 
    E.g., for a marker named by SOX10, SOX10_positive will be added, which has two values: True or False.

    @adata: AnnData
    @marker: target_marker
    @gate: the gate value for the given marker
    '''
    
    newkey = target_marker+"_positive"
    adata.obs[newkey] = (adata[list(adata.obs_names), [target_marker]].X >= gate).flatten().tolist()
    if debug: print(adata.obs[newkey].value_counts())

