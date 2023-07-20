import shapely.geometry as geometry
from descartes import PolygonPatch
import matplotlib.pyplot as plt
import alphashape
import scipy as sp
import networkx as nx
import numpy as np

def compute_roundness(polygon):
    """[summary]

    Args:
        polygon (shapely.geometry.polygon object): polygon of a cluster

    Returns:
        float: roundness of the target cluster
    """
    alpha_shape = polygon
    roundness = 4*sp.pi*alpha_shape.area/pow(alpha_shape.length,2)
    return roundness
