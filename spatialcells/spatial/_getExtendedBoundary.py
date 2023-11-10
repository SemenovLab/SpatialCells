from ._bufferBoundary import bufferBoundary


def getExtendedBoundary(boundaries, offset=200):
    """
    Get the expand the boundary of the communities of interest
    
    :param boundaries: the boundaries of components
    :param offset: the offset to extend the boundary
    :returns: the extended boundaries of components.
    """
    return bufferBoundary(boundaries, offset)
