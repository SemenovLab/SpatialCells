from ._utils import *


def getShrunkenBoundary(boundaries, offset=200, minsize=30):
    """
    Get the shrunken boundary of the communities of interest
    :param boundaries: the boundaries of components
    :param offset: the offset to shrink the boundary
    :param minsize: Components with size less than minsize will be removed
    :return: the shrunken boundaries of components.
    """
    return getBufferedBoundary(boundaries, -offset, minsize)
