import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from ._utils import *

def getExtendedBoundary(boundaries, offset=200, minsize=30):
    return getBufferedBoundary(boundaries, offset, minsize)





