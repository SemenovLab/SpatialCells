

def calculate_area(polygon_in):
    """[summary]

    Args:
        polygon_in (shapely.geometry.polygon object): polygon of a cluster

    Returns:
       float: area for polygon of a cluster
    """
    x, y = zip(*list(polygon_in.exterior.coords))
    return abs(sum(x[i-1]*y[i]-x[i]*y[i-1] for i in range(len(x)))) / 2.