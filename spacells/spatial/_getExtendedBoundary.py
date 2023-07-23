import numpy as np
import math
from ._utils import *

def getExtendedBoundary(boundary, offset=400, alpha=200):
    '''
    '''
    # print("get points on edges")
    boundary_points = boundary[:,0,:]

    # print("get stretched points")
    pointsStretched = bufferPoints(boundary_points, offset, n=50)
    
    # print("get edge_components of stretched points")
    # set of (i,j) point pairs representing edges of the alpha-shape.
    edge_components = getAlphaShapes(pointsStretched, alpha)
    # print("edge_components:", len(edge_components), 
    #     [len(comp) for comp in edge_components])
    
    extended_boundaries = []
    for comp in edge_components:
        keep = True
        # print(len(comp))
        for i in range(len(comp)):
            if isInRegion(comp[i][0], boundary):
                # print(i, "in", comp[i][0])
                keep = False
                break
        if keep:
            extended_boundaries.append(comp)
    # if len(extended_boundaries) == 0:
    #     print("wrong?")
    return np.concatenate(extended_boundaries)

