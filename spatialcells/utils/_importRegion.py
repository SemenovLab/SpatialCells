import pandas as pd
from shapely import from_geojson, make_valid
from shapely.geometry import Polygon, MultiPolygon


def importRegion(filename, scale=1.0):
    """Import a region from a file. The file type is determined by the 
    file extension.

    :param filename: The name of the file to import from.
    :param scale: The scale factor to apply to the coordinates.
    :return: region/boundary object
    """
    filetype = filename.split(".")[-1]
    if filetype in ["geojson", "json"]:
        return _from_geojson(filename)
    elif filetype == "csv":
        return _from_csv(filename, scale)
    else:
        raise ValueError("File type not supported. Please use geojson or csv.")
    
def _from_csv(filename, scale):
    df = pd.read_csv(filename)
    polygons = []
    poly_coords = df["all_points"].tolist()
    for coords in poly_coords:
        coords = list(map(
            lambda x:tuple(map(
                lambda y: float(y) * scale, 
                x.split(",")
            )), 
            coords.split()
        ))
        polygons.append(Polygon(coords))
    return make_valid(MultiPolygon(polygons))

def _from_geojson(filename):
    with open(filename, "r") as f:
        geojson = f.read()
    return from_geojson(geojson)
