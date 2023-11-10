import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects

from ._plotBoundary import plotBoundary


def plotRegions(
    regions_list,
    ax=None,
    add_label=True,
    label_prefix="",
    x_offset=0,
    y_offset=0,
    fontsize=12,
    label_bounds=None,
    colors_list=["k", "r", "b", "g"],
    fontcolor="black",
    foreground_linewidth=5,
    **kwargs
):
    """
    Plot the regions in the regions_list and label them with their index in
    the list.
    labels are placed at the lower left corner of each region. Their position
    can be adjusted by x_offset and y_offset, and the overall label_bounds.
    
    :param regions_list: a list of MultiPolygon regions
    :param add_label: whether to label the regions
    :param ax: the matplotlib ax object.
    :param label_prefix: a prefix to the label of each region
    :param x_offset: the x offset of each label
    :param y_offset: the y offset of the label
    :param fontsize: the fontsize of the label
    :param label_bounds: the bounds of all the label (minx, miny, maxx, maxy).
        If None, the bounds of the regions_list will be used.
    :param colors_list: a list of colors that will be cycled through for each
        region.
    :param fontcolor: the fontcolor of the label
    :param foreground_linewidth: the foreground linewidth of the label. Adjusts
        the thickness of the white border around each label.
    :param kwargs: kwargs for matplotlib.pyplot.plot
    """
    if ax is None:
        ax = plt
    for i, region in enumerate(regions_list):
        bounds = list(region.bounds)
        if label_bounds is not None:
            bounds[0] = max(bounds[0], label_bounds[0])
            bounds[1] = max(bounds[1], label_bounds[1])
            bounds[2] = min(bounds[2], label_bounds[2])
            bounds[3] = min(bounds[3], label_bounds[3])
        plotBoundary(region, ax=ax, color=colors_list[i % len(colors_list)], **kwargs)
        if not add_label:
            continue
        label = label_prefix + str(i)
        txt = ax.text(
            bounds[0] + x_offset,
            bounds[1] + y_offset,
            label,
            fontsize=fontsize,
            color=fontcolor,
        )
        txt.set_path_effects(
            [PathEffects.withStroke(linewidth=foreground_linewidth, foreground="w")]
        )
