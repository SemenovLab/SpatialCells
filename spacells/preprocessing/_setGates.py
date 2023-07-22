import numpy as np
import pandas as pd

def setGates(adata, gate_file, debug=True):
    '''

    @adata: AnnData
    @gate_file: 
    '''
    
    manual_gate = pd.read_csv(gate_file)
    
    # Some markers may not be in adata
    markers_gated = []
    for m in list(adata.var_names):
        if m in list(manual_gate["markers"]):
            markers_gated.append(m)
    
    # Set gates
    for target_marker in markers_gated:
        val = manual_gate.loc[manual_gate['markers'] == target_marker, 'gate'].values[0]
        newkey = target_marker+"_b"
        adata.obs[newkey] = (adata[list(adata.obs_names), [target_marker]].X >= val).flatten().tolist()
        if debug: print("gate:", val, adata.obs[newkey].value_counts())

