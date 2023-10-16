def getRegionCentroid(boundary):
    """
    Get the centroid of a region defined by a list of region boundary components.
    :param boundary: A MultiPolygon object defining the boundary of the region
    :return: The centroid of the region
    """
    xy = boundary.centroid.xy
    return [xy[0][0], xy[1][0]]
