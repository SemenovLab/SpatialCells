import numpy as np
import math
from _utils import *

def getShrunkenBoundary(boundary, offset=200, alpha=100, minsize=30):

    # print("get points on edges")
    boundary_points = boundary[:,0,:]

    # print("get stretched points")
    pointsStretched = bufferPoints(boundary_points, offset, n=100)

    # print("get edge_components of stretched points")
    edge_components = getAlphaShapes(pointsStretched, alpha)
    # print("edge_components:", len(edge_components), 
    #     [len(comp) for comp in edge_components])

    shrunken_boundaries = []
    for comp in edge_components:
        keep = True
        # print(len(comp))
        for i in range(len(comp)):
            if len(comp) < minsize:
                keep = False
                break
            point = comp[i][0]
            if not isInRegion(point, boundary):
                # print("out", point)
                keep = False
                break
        if keep:
            shrunken_boundaries.append(comp)
    return np.concatenate(shrunken_boundaries)


