from ._utils import *


def getExtendedBoundary(boundaries, offset=200, minsize=30):
    """
    Get the extended boundary of the communities of interest
    :param boundaries: the boundaries of components
    :param offset: the offset to extend the boundary
    :param minsize: Components with size less than minsize will be removed
    :return: the extended boundaries of components.
    """
    return getBufferedBoundary(boundaries, offset, minsize)
