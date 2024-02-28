import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import numpy as np


def plotCommunities(adata, ret, communitycolumn, plot_first_n_clusters=10, **kwargs):
    """
    Plot largest communities on a scatter plot. Each community is labelled as:
    (rank in descending number of cells : index of community)

    :param adata: AnnData object
    :param ret: return value of spc.spa.getCommunities
    :param communitycolumn: column name of the community column in adata.obs
    :param plot_first_n_clusters: plot and number the largest n communities
    :param kwargs: keyword arguments for matplotlib.pyplot.plot
    """
    if "s" not in kwargs:
        kwargs["s"] = 2
    markersize = kwargs["s"]
    kwargs.pop("s", None)
    if "fontsize" not in kwargs:
        kwargs["fontsize"] = 12
    fontsize = kwargs["fontsize"]
    kwargs.pop("fontsize", None)
    if "ax" not in kwargs:
        fig, ax = plt.subplots(figsize=(10, 8))
    else:
        ax = kwargs["ax"]
        kwargs.pop("ax", None)

    labels_sorted, db = ret

    adata_tmp = adata[adata.obs[communitycolumn] != -2]
    X = adata_tmp.obs[["X_centroid", "Y_centroid"]].to_numpy()
    labels = db.labels_
    unique_labels = set(labels)
    core_samples_mask = np.zeros_like(labels, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]

    ## all points
    ax.scatter(
        *zip(*adata.obs[["X_centroid", "Y_centroid"]].to_numpy()),
        s=3,
        color="grey",
        alpha=0.2
    )

    # Points of interest: outliers
    ax.scatter(X[:, 0], X[:, 1], alpha=0.5, color="orange", s=markersize / 2)

    # clusters
    idx = 0
    for npoints, k in labels_sorted:
        if k == -1:
            continue
        if idx >= plot_first_n_clusters:
            break

        col = colors[idx]
        class_member_mask = labels == k
        mask = class_member_mask & core_samples_mask
        xy = X[mask]
        plt.plot(
            xy[:, 0],
            xy[:, 1],
            "o",
            markerfacecolor=tuple(col),
            markeredgecolor=tuple(col),
            markersize=markersize,
            **kwargs
        )
        # plt.annotate(str(idx + 1) + ":" + str(k), (xy[0, 0], xy[0, 1]))
        txt = plt.text(
            xy[0, 0],
            xy[0, 1],
            str(idx + 1) + ":" + str(k),
            fontsize=fontsize,
            color="black",
        )
        txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground="w")])
        idx += 1
