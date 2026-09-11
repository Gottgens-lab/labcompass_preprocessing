"""Figure helpers used by the preprocessing notebooks."""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb

__all__ = [
    "plot_marker",
    "plot_leiden_clusters",
    "plot_ab_histograms",
    "plot_umap_by_group",
]


def plot_marker(
    adata,
    markers,
    threshold=2,
    n_cols=6,
    layer=None,
):
    """Draw one UMAP panel per marker, colouring cells above `threshold` in red.

    The intensity is read from ``adata.raw`` when `layer` is None, which is where
    the notebooks keep the background-subtracted, zero-to-one scaled values.
    Pass a layer name to read from ``adata.layers`` instead.

    `threshold` is in the units of whichever matrix is read.  For the
    zero-to-one scaled "positives" values the notebooks use 0.01.
    """
    if isinstance(markers, str):
        markers = [markers]

    if layer is None:
        if adata.raw is None:
            raise ValueError(
                "adata.raw is not set; either set it or pass layer=..."
            )
        source_var_names = adata.raw.var_names
        source_X = adata.raw.X
    else:
        source_var_names = adata.var_names
        source_X = adata.layers[layer]

    missing = [m for m in markers if m not in source_var_names]
    if missing:
        raise KeyError(f"markers not present in the intensity matrix: {missing}")

    n_rows = -(len(markers) // -n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3 * n_cols, 3 * n_rows))
    axes = np.atleast_1d(axes).flatten()

    umap = adata.obsm["X_umap"]

    for i, marker in enumerate(markers):
        marker_idx = source_var_names.get_loc(marker)
        values = np.asarray(source_X[:, marker_idx]).ravel()
        mask = values > threshold

        axes[i].scatter(
            umap[~mask, 0], umap[~mask, 1], color="lightgrey", alpha=0.5, s=0.02
        )
        axes[i].scatter(umap[mask, 0], umap[mask, 1], color="red", alpha=1, s=0.1)
        axes[i].set_title(f"Marker: {marker}")
        axes[i].set_xlabel("UMAP1")
        axes[i].set_ylabel("UMAP2")
        axes[i].set_box_aspect(1)

    for j in range(len(markers), len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()
    return fig


def plot_leiden_clusters(adata, cluster_key="leiden", n_cols=4, palette=None):
    """Draw one UMAP panel per cluster, colouring that cluster's cells.

    `palette` maps a cluster label to a colour.  When it is None the clusters
    are coloured by the ``tab20`` colormap, which repeats after twenty clusters.
    """
    clusters = adata.obs[cluster_key].astype(str)
    unique_clusters = list(clusters.unique())

    if palette is None:
        cmap = plt.get_cmap("tab20")
        palette = {cl: cmap(i % 20) for i, cl in enumerate(unique_clusters)}

    n_rows = -(len(unique_clusters) // -n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3 * n_cols, 3 * n_rows))
    axes = np.atleast_1d(axes).flatten()

    umap = adata.obsm["X_umap"]

    for i, cluster in enumerate(unique_clusters):
        mask = (clusters == cluster).values
        axes[i].scatter(
            umap[~mask, 0], umap[~mask, 1], color="lightgrey", alpha=0.5, s=0.02
        )
        axes[i].scatter(
            umap[mask, 0], umap[mask, 1], color=palette[cluster], alpha=1, s=0.1
        )
        axes[i].set_title(f"Cluster: {cluster}")
        axes[i].set_xlabel("UMAP1")
        axes[i].set_ylabel("UMAP2")
        axes[i].set_box_aspect(1)

    for j in range(len(unique_clusters), len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()
    return fig


def plot_ab_histograms(adata, ab_column="antibody", bins=40, n_cols=4):
    """Draw one intensity histogram per antibody of `adata`.

    `ab_column` names the column of ``adata.var`` that holds the antibody name.
    """
    antibodies = adata.var[ab_column].astype(str).tolist()

    n_rows = -(len(antibodies) // -n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3 * n_rows))
    axes = np.atleast_1d(axes).flatten()

    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()

    for i, antibody in enumerate(antibodies):
        values = np.asarray(X[:, i]).ravel()
        sb.histplot(values, bins=bins, ax=axes[i])
        axes[i].set_title(antibody)
        axes[i].set_xlabel("intensity")

    for j in range(len(antibodies), len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()
    return fig


def plot_umap_by_group(adata, obs="experiment_number", n_cols=3, size=5, cmap="tab10"):
    """Draw one UMAP panel per value of ``adata.obs[obs]``.

    Every panel shows all the cells in light grey and the cells of that one group
    in colour, so a group that sits in only part of the embedding is visible at a
    glance.  Groups are sorted, so the panel order does not depend on the order
    in which the wells were read.
    """
    groups = sorted(adata.obs[obs].unique().tolist())
    n_rows = -(len(groups) // -n_cols)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
    axes = np.atleast_1d(axes).flatten()

    umap = adata.obsm["X_umap"]
    colours = plt.get_cmap(cmap)(np.linspace(0, 1, max(len(groups), 2)))

    for i, group in enumerate(groups):
        mask = (adata.obs[obs] == group).values
        axes[i].scatter(umap[:, 0], umap[:, 1], s=size, color="lightgrey", linewidths=0)
        axes[i].scatter(
            umap[mask, 0], umap[mask, 1], s=size, color=colours[i], linewidths=0
        )
        axes[i].set_title(str(group))
        axes[i].set_xticks([])
        axes[i].set_yticks([])

    for j in range(len(groups), len(axes)):
        axes[j].axis("off")

    plt.tight_layout()
    return fig
