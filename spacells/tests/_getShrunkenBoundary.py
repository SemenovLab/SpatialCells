import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from _utils import *

def getExtendedBoundary(boundaries, offset=200, minsize=30):
    return getBufferedBoundary(boundaries, offset, minsize)

def getShrunkenBoundary(boundaries, offset=200, minsize=30):
    return getBufferedBoundary(boundaries, -offset, minsize)

def getBufferedBoundary(boundaries, offset=200, minsize=30):
    buffered_boundaries = []
    outer_multipolygon = Polygon()
    inner_multipolygon = Polygon()
    for boundary_set in boundaries:
        for i, boundary in enumerate(boundary_set):
            boundary_points = boundary[:, 0, :]
            boundary_polygon = Polygon(boundary_points)
            if i == 0:
                outer_multipolygon = outer_multipolygon | boundary_polygon.buffer(-offset)
            else:
                inner_multipolygon = inner_multipolygon | boundary_polygon.buffer(offset)
    s_boundary_polygons = outer_multipolygon - inner_multipolygon
    if isinstance(s_boundary_polygons, MultiPolygon):
        s_boundary_polygons = s_boundary_polygons.geoms
    else:
        s_boundary_polygons = [s_boundary_polygons]
    for s_boundary_polygon in s_boundary_polygons:
        for ring in [s_boundary_polygon.exterior] + list(s_boundary_polygon.interiors):
            buffered_points = np.array(ring.coords)
            if buffered_points.shape[0] < minsize:
                continue
            buffered_edges = np.stack([np.roll(buffered_points, -1, axis=0), buffered_points], axis=1)
            buffered_boundaries.append(buffered_edges)
    return [buffered_boundaries]