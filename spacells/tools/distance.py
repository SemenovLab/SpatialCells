import networkx as nx
import numpy as np

def compute_farthest_two_points(polygon):
    """[summary]

    Args:

        polygon (shapely.geometry.polygon object): polygon of a cluster

    Returns:
        tuple(float, (float, float)): (maximum distance, coordinates of two farthest points in th cluster)
    """
    alpha_shape = polygon
    G = nx.complete_graph(set(alpha_shape.exterior.coords))
    max_result = ()
    max_dist = float('-inf')
    for i in G.edges():
        tmp = np.linalg.norm(np.asarray(i[0]) - np.asarray(i[1]))
        if tmp > max_dist:
            max_dist = tmp
            max_result = (i[0], i[1])
    return (max_dist, max_result)