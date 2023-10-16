from ._bufferBoundary import bufferBoundary


def getShrunkenBoundary(boundaries, offset=200):
    """
    Shrink the boundary of the communities of interest.
    Wrapper of bufferBoundary with negative offset.
    :param boundaries: the boundaries of components
    :param offset: the offset to shrink the boundary
    :return: the shrunken boundaries of components.
    """
    return bufferBoundary(boundaries, -offset)
