import pandas as pd
from shapely import to_geojson


def exportRegion(region, filename):
    """Export a region to a file. The file type is determined by the
    file extension.

    :param region: The region/boundary object to export.
    :param filename: The name of the file to export to.
    :return: geojson string
    """
    filetype = filename.split(".")[-1]
    if filetype in ["geojson", "json"]:
        return _to_geojson(region, filename)
    elif filetype == "csv":
        return _to_csv(region, filename)
    else:
        raise ValueError("File type not supported. Please use geojson or csv.")
    
def _to_csv(region, filename):
    poly_coords = []
    for poly in region.geoms:
        poly = poly.exterior
        coords = " ".join(map(lambda x: ",".join(map(str, x)), poly.coords))
        poly_coords.append(coords)
    df = pd.DataFrame({
        "Id": range(len(poly_coords)),
        "Name": ["Polygon"] * len(poly_coords),
        "Text": [""] * len(poly_coords),
        "type": ["Polygon"] * len(poly_coords),
        "all_points": poly_coords,
        "X": ["-1"] * len(poly_coords),
        "Y": ["-1"] * len(poly_coords),
        "RadiusX": ["-1"] * len(poly_coords),
        "RadiusY": ["-1"] * len(poly_coords),
        "Width": ["-1"] * len(poly_coords),
        "Height": ["-1"] * len(poly_coords),
        "all_transforms":[ "-1"] * len(poly_coords),
    })
    df.to_csv(filename, index=False)
    return df
        
    
def _to_geojson(region, filename):
    geojson = to_geojson(region)
    with open(filename, "w") as f:
        f.write(geojson)
    return geojson
