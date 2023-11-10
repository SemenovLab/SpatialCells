from shapely.geometry import MultiPolygon
from shapely.validation import make_valid


def bufferBoundary(boundary, offset):
    """
    Buffer a boundary by a given offset. Negative offset will shrink the boundary.
    
    :param boundary: the boundary to be buffered
    :param offset: the offset
    :returns: the buffered boundary
    """
    new_boundary = boundary.buffer(offset)
    if new_boundary.geom_type == "Polygon":
        new_boundary = MultiPolygon([new_boundary])
    return make_valid(new_boundary)
