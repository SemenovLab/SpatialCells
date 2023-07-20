
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.sparse import csr_matrix

import networkx as nx

def compute_max_distance(csr_matrix_input):
    """[summary]

    Args:
        csr_matrix_input (scipy.sparse.csr_matrix object): graph of the target cluster

    Returns:
        int: maximum distance calculated from MST with unit length, 1 for each edge
    """
    MST_result = minimum_spanning_tree(csr_matrix_input)
    Max_length = csr_matrix.sum(MST_result)
    return Max_length

def compute_mst_distance(graph):
    """[summary]
    Args:
        graph: graph of the target cluster
    Returns:
        int: maximum distance calculated from MST with unit length, 1 for each edge
    """
    T = nx.minimum_spanning_tree(graph)
    return len(T.edges)
    