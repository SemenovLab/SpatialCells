import shapely
from shapely.geometry import Polygon


def getRegionArea(boundary, exclude_holes=True):
    """
    Get the area of a region defined by a MultiPolygon boundary.
    If exclude_holes is True, the area of the holes in the region is
    subtracted from the overall area of the region.
    :param boundary: MultiPolygon boundary of the region
    :param exclude_holes: whether to exclude the holes in the region
    :return: area of the region
    """
    if exclude_holes:
        return boundary.area
    else:
        extern_polygons = [Polygon(polygon.exterior) for polygon in boundary.geoms]
        overall_polygon = shapely.unary_union(extern_polygons)
        return overall_polygon.area
