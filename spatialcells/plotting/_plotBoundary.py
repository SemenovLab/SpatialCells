import matplotlib.pyplot as plt


def plotBoundary(boundary, ax=None, **kwargs):
    """
    Plot the boundary
    :param boundary: the MultiPolygon boundary
    :param ax: the matplotlib ax object.
    :param kwargs: kwargs for matplotlib.pyplot.plot
    """
    if ax is None:
        ax = plt
    if "c" not in kwargs and "color" not in kwargs:
        kwargs["c"] = "k"
    for polygon in boundary.geoms:
        ax.plot(*polygon.exterior.xy, **kwargs)
        kwargs.pop("label", None)
        for interior in polygon.interiors:
            ax.plot(*interior.xy, **kwargs)
