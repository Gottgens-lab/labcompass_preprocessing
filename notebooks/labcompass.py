"""Shared helpers for the LabCompass preprocessing notebooks.

Only the long, genuinely shared pieces live here.  Everything short enough to
read at a glance is written out in the notebook cell that uses it, so a reader
can follow a stage without opening this file.

The notebooks import it by sitting one directory below it::

    import sys
    sys.path.append("..")
    import labcompass as lc

What is here, and how many of the 32 notebooks import each piece:

| function | notebooks | why it is not written out in the cell |
|---|---|---|
| `downsample` | 27 | repeated in almost every notebook |
| `apply_saved_gates` | 18 | repeated, and it reads a file format |
| `get_axis_data`, `read_gate_file` | 18 | `apply_saved_gates` and `DrawGate` both need them |
| `clean_intensities` | 9 | 38 lines, and the order of its steps decides the result |
| `logicle_transform` | 9 | 24 lines, and it checks the parameter lengths |
| `compute_umap` | 9 | 21 lines of scanpy calls on a scaled copy |
| `plot_marker`, `plot_umap_by_group` | 9 | 43 and 24 lines of plotting |
| `DrawGate` | 2 | 234 lines of interactive matplotlib |
| `rotate_umap`, `plot_leiden_clusters` | 2 | repeated in both loops that cluster |
"""

from collections import Counter, defaultdict

import anndata as ann
import matplotlib.pyplot as plt
import numpy as np
import pytometry as pm
import scanpy as sc
import seaborn as sb
from matplotlib.patches import Polygon
from matplotlib.path import Path
from scipy import sparse
from scipy.spatial import ConvexHull, Delaunay

__all__ = [
    "get_axis_data",
    "read_gate_file",
    "apply_saved_gates",
    "DrawGate",
    "downsample",
    "clean_intensities",
    "logicle_transform",
    "compute_umap",
    "rotate_umap",
    "plot_marker",
    "plot_leiden_clusters",
    "plot_ab_histograms",
    "plot_umap_by_group",
]


# =====================================================================
# Gating
# =====================================================================
#
# A gate file is a ``.npz`` archive with three parallel arrays.  Entry ``i`` of
# ``polygons`` is an ``(n_vertices, 2)`` array of vertex coordinates, and entries
# ``i`` of ``x_cols`` and ``y_cols`` name the two axes the polygon was drawn
# against.  An axis name is looked up first in ``adata.obs`` and then in
# ``adata.var_names``, so a gate may be drawn against a scatter channel such as
# ``"FSC-A :: FSC - Area"``, a derived observation such as ``"7-AAD"``, or a
# marker channel such as ``"CD235a : AF700 - Area"``.


def get_axis_data(adata, key):
    """Return a one-dimensional array for `key`, from ``.obs`` or from ``.X``.

    `key` is looked up in ``adata.obs.columns`` first and in ``adata.var_names``
    second.  A `KeyError` is raised when it is in neither, which is what makes a
    typo in a gate file fail loudly instead of producing an empty gate.
    """
    if key in adata.obs.columns:
        return adata.obs[key].values

    if key in adata.var_names:
        values = adata[:, key].X
        if hasattr(values, "toarray"):
            values = values.toarray()
        return np.asarray(values).ravel()

    raise KeyError(f"{key!r} not found in adata.obs or adata.var_names")


def read_gate_file(gate_file):
    """Return `(polygons, x_cols, y_cols)` from a gate file, as plain lists.

    Each polygon comes back as a float array of shape ``(n_vertices, 2)``.
    `np.savez` stores a list of equal-length polygons as one regular
    three-dimensional array rather than as an array of objects, so the polygons
    are rebuilt element by element to cover both storage shapes.
    """
    data = np.load(gate_file, allow_pickle=True)
    polygons = [np.asarray(p, dtype=float) for p in data["polygons"]]
    x_cols = [str(c) for c in data["x_cols"]]
    y_cols = [str(c) for c in data["y_cols"]]

    if not (len(polygons) == len(x_cols) == len(y_cols)):
        raise ValueError(
            f"Mismatch in {gate_file}: {len(polygons)} polygons, "
            f"{len(x_cols)} x-columns, {len(y_cols)} y-columns"
        )
    return polygons, x_cols, y_cols


def apply_saved_gates(adata, gate_file, mode="all"):
    """Apply every polygon in `gate_file` to `adata` and return the kept cells.

    With ``mode="all"`` a cell is kept only when it falls inside every polygon,
    which is the sequential gating strategy the notebooks use.  With
    ``mode="any"`` a cell is kept when it falls inside at least one polygon.

    Returns the gated copy of `adata` and the list of per-gate boolean masks,
    each mask having one entry per cell of the input.
    """
    polygons, x_cols, y_cols = read_gate_file(gate_file)

    if not polygons:
        raise ValueError(f"{gate_file} contains no polygons")

    masks = []
    for poly, x, y in zip(polygons, x_cols, y_cols):
        pts = np.c_[get_axis_data(adata, x), get_axis_data(adata, y)]
        masks.append(Path(poly).contains_points(pts))

    if mode == "any":
        final_mask = np.logical_or.reduce(masks)
    elif mode == "all":
        final_mask = np.logical_and.reduce(masks)
    else:
        raise ValueError("mode must be 'any' or 'all'")

    return adata[final_mask].copy(), masks


class DrawGate:
    """An interactive polygon gate over two channels of an AnnData object.

    Drag the red vertex handles to reshape the polygon, close the figure, then
    call `filter_cells` to subset and `save_gate` to append the polygon to a
    gate file.  `x` and `y` follow the lookup rule of `get_axis_data`.

    Needs an interactive matplotlib backend, so it only works inside Jupyter
    with ``%matplotlib widget`` and `ipympl` installed.

    Set `high_res` to plot every event with `scanpy.pl.scatter`, which is exact
    but slow above a few hundred thousand cells.  Leave it False to plot a
    density image with `pytometry.pl.scatter_density` instead.

    `auto_hull` may be False for a regular polygon of `n_vertices` vertices
    placed in the middle of the axes, ``"rect"`` for the axis-aligned bounding
    rectangle of the data, or ``"convex"`` for an alpha shape around the data.
    """

    def __init__(
        self,
        adata,
        x="FSC-A :: FSC - Area",
        y="SSC-A :: SSC - Area",
        n_vertices=5,
        xlim=None,
        ylim=None,
        xscale="linear",
        yscale="linear",
        auto_hull=False,
        high_res=False,
    ):
        self.adata = adata
        self.x = x
        self.y = y

        self.fig, self.ax = plt.subplots(figsize=(6, 6))

        if high_res:
            sc.pl.scatter(adata, x=x, y=y, ax=self.ax, color="blue", show=False)
            if xlim is not None:
                self.ax.set_xlim(xlim)
            if ylim is not None:
                self.ax.set_ylim(ylim)
        else:
            pm.pl.scatter_density(
                adata,
                x=x,
                y=y,
                x_lim=xlim,
                y_lim=ylim,
                ax=self.ax,
                y_scale=yscale,
                x_scale=xscale,
            )

        self.polys = []
        self.verts = []
        self.handles = []
        self.dragging = None

        pts = np.c_[get_axis_data(adata, x), get_axis_data(adata, y)]

        if auto_hull == "convex":
            poly = self._alpha_shape_polygon(pts, alpha=1e-6)
        elif auto_hull == "rect":
            poly = self._bounding_rect(pts)
        else:
            xmin, xmax = self.ax.get_xlim()
            ymin, ymax = self.ax.get_ylim()
            cx = (xmin + xmax) / 2
            cy = (ymin + ymax) / 2
            scale = 0.15 * min(xmax - xmin, ymax - ymin)
            theta = np.linspace(0, 2 * np.pi, n_vertices, endpoint=False)
            poly = np.c_[cx + scale * np.cos(theta), cy + scale * np.sin(theta)]

        patch = Polygon(poly, closed=True, fill=False, edgecolor="red", lw=2)
        self.ax.add_patch(patch)
        handle = self.ax.plot(poly[:, 0], poly[:, 1], "ro", ms=6)[0]

        self.polys.append(patch)
        self.verts.append(poly)
        self.handles.append(handle)

        self.fig.canvas.mpl_connect("button_press_event", self.on_press)
        self.fig.canvas.mpl_connect("motion_notify_event", self.on_motion)
        self.fig.canvas.mpl_connect("button_release_event", self.on_release)

        self.ax.set_title("Drag the red dots to reshape the gate, then close")

    # ---------------- mouse interaction ----------------

    def _pixel_distance(self, x, y, xv, yv):
        xy_disp = self.ax.transData.transform(np.c_[xv, yv])
        x_disp, y_disp = self.ax.transData.transform((x, y))
        return np.hypot(xy_disp[:, 0] - x_disp, xy_disp[:, 1] - y_disp)

    def on_press(self, event):
        if event.inaxes != self.ax or event.xdata is None:
            return
        for gi, poly in enumerate(self.verts):
            d = self._pixel_distance(event.xdata, event.ydata, poly[:, 0], poly[:, 1])
            vi = int(np.argmin(d))
            if d[vi] < 10:
                self.dragging = (gi, vi)
                return

    def on_motion(self, event):
        if self.dragging is None or event.inaxes != self.ax:
            return
        if event.xdata is None or event.ydata is None:
            return
        gi, vi = self.dragging
        self.verts[gi][vi] = [event.xdata, event.ydata]
        self.polys[gi].set_xy(self.verts[gi])
        self.handles[gi].set_data(self.verts[gi][:, 0], self.verts[gi][:, 1])
        self.fig.canvas.draw_idle()

    def on_release(self, event):
        self.dragging = None

    # ---------------- results ----------------

    def get_masks(self):
        """Return one boolean mask per polygon, over the cells of `self.adata`."""
        pts = np.c_[
            get_axis_data(self.adata, self.x), get_axis_data(self.adata, self.y)
        ]
        return [Path(poly).contains_points(pts) for poly in self.verts]

    def filter_cells(self, mode="any"):
        """Return the cells this gate keeps, as a new AnnData object."""
        masks = self.get_masks()
        if mode == "any":
            final = np.logical_or.reduce(masks)
        else:
            final = np.logical_and.reduce(masks)
        return self.adata[final].copy()

    def save_gate(self, filename):
        """Append this gate's polygon and its two axis names to `filename`.

        The gate file is created when it does not exist, so starting a gating
        strategy with an empty gate file is what keeps a rerun from stacking new
        gates on top of the gates of the previous run.
        """
        import os

        if os.path.exists(filename):
            polygons, xs, ys = read_gate_file(filename)
        else:
            polygons, xs, ys = [], [], []

        polygons.append(np.asarray(self.verts[0], dtype=float))
        xs.append(self.x)
        ys.append(self.y)

        np.savez(
            filename,
            polygons=np.array(polygons, dtype=object),
            x_cols=np.array(xs, dtype=object),
            y_cols=np.array(ys, dtype=object),
        )

    # ---------------- automatic starting polygons ----------------

    @staticmethod
    def _bounding_rect(pts):
        """Return the four vertices of the axis-aligned bounding box of `pts`."""
        pts = np.asarray(pts)
        xmin, ymin = pts.min(axis=0)
        xmax, ymax = pts.max(axis=0)
        return np.array([[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]])

    @staticmethod
    def _alpha_shape_polygon(points, alpha):
        """Return an alpha shape around `points`, falling back to a convex hull.

        A triangle of the Delaunay triangulation is kept when its circumradius is
        below ``1 / alpha``.  The boundary edges of the kept triangles are then
        walked to give one closed ring.  The convex hull is returned instead
        whenever no triangle is kept, fewer than three boundary edges survive, or
        the walk cannot be completed.
        """
        pts = np.asarray(points)
        tri = Delaunay(pts)

        edges = []
        for simplex in tri.simplices:
            pa, pb, pc = pts[simplex]
            a = np.linalg.norm(pb - pa)
            b = np.linalg.norm(pc - pb)
            c = np.linalg.norm(pa - pc)
            s = (a + b + c) / 2
            area2 = s * (s - a) * (s - b) * (s - c)
            if area2 <= 0:
                continue
            radius = a * b * c / (4.0 * np.sqrt(area2))
            if radius < 1.0 / alpha:
                edges.extend(
                    [
                        tuple(sorted((simplex[0], simplex[1]))),
                        tuple(sorted((simplex[1], simplex[2]))),
                        tuple(sorted((simplex[2], simplex[0]))),
                    ]
                )

        if not edges:
            return pts[ConvexHull(pts).vertices]

        counts = Counter(edges)
        boundary = [e for e in counts if counts[e] == 1]
        if len(boundary) < 3:
            return pts[ConvexHull(pts).vertices]

        adj = defaultdict(list)
        for i, j in boundary:
            adj[i].append(j)
            adj[j].append(i)

        start = next((k for k, v in adj.items() if len(v) == 1), boundary[0][0])

        ring = [start]
        visited = {start}
        prev = None
        current = start
        # The walk is bounded by the number of boundary vertices, so a branching
        # or broken boundary ends the walk instead of looping for ever.
        for _ in range(len(adj)):
            neighbours = [n for n in adj[current] if n != prev]
            if not neighbours:
                break
            nxt = neighbours[0]
            if nxt == start or nxt in visited:
                break
            ring.append(nxt)
            visited.add(nxt)
            prev, current = current, nxt

        if len(ring) < 3:
            return pts[ConvexHull(pts).vertices]

        return pts[np.array(ring, dtype=int)]


# =====================================================================
# Sampling
# =====================================================================


def downsample(adata, obs="experiment_number", target_size=None, seed=0):
    """Cap every group in ``adata.obs[obs]`` at `target_size` cells.

    A group with more than `target_size` cells is sampled without replacement.
    A group with fewer cells is kept whole, so the groups are only balanced when
    every group is at least `target_size` cells large.  When `target_size` is
    None the size of the smallest group is used, which does balance the groups.

    Returns a new AnnData object.  The input is not modified.

    Reproducibility: this calls ``np.random.seed(seed)``, which reseeds the
    global NumPy generator as a side effect.  Two calls with the same seed, the
    same group sizes and the same NumPy version select the same cells.  Any code
    that draws random numbers after this call is affected by the reseeding.
    """
    np.random.seed(seed)

    group_sizes = adata.obs[obs].value_counts()

    if target_size is None:
        target_size = int(group_sizes.min())

    keep_positions = []
    obs_values = adata.obs[obs].values

    for group in group_sizes.index:
        group_positions = np.flatnonzero(obs_values == group)
        if len(group_positions) > target_size:
            chosen = np.random.choice(
                len(group_positions), size=target_size, replace=False
            )
            keep_positions.append(group_positions[chosen])
        else:
            keep_positions.append(group_positions)

    return adata[np.concatenate(keep_positions)].copy()


# =====================================================================
# Cleaning and transformation
# =====================================================================


def clean_intensities(
    adata,
    marker_cutoff,
    percentile=False,
    p=99.9997,
    subtract_fmo_cutoff=False,
    zero_qc=False,
    add_fmo_cutoff=False,
    negatives_to_zeroes=False,
    shift=False,
    scaling=False,
):
    """Apply the optional cleaning steps to ``adata.X``, in place.

    The steps run in this fixed order, and each one is skipped unless its flag
    is set.  The order is the point of this function: a different order gives a
    different result.

    1. `percentile`: drop a cell when any channel is above the `p`-th percentile
       of that channel.  This removes the extreme upper tail.  Cells are
       removed, so ``adata.obs`` shrinks with ``adata.X``.
    2. `subtract_fmo_cutoff`: subtract the per-antibody cutoff from every
       intensity, which puts the fluorescence-minus-one background at zero.
    3. `zero_qc`: drop a cell when no channel is above zero.  Run after the
       subtraction, this drops cells with no signal in any marker.
    4. `add_fmo_cutoff`: add the per-antibody cutoff back on.
    5. `negatives_to_zeroes`: replace every negative intensity with zero.
    6. `shift`: subtract the per-channel minimum.
    7. `scaling`: divide each channel by its maximum, putting it on a zero to
       one scale.  A channel whose maximum is zero is left alone.

    `marker_cutoff` maps an antibody name to its cutoff, and is only read when
    `subtract_fmo_cutoff` or `add_fmo_cutoff` is set.  Every name in
    ``adata.var_names`` must be a key of it, or a `KeyError` is raised.

    Returns None.  `adata` is modified in place.
    """
    # The dtype is deliberately not forced.  The intensities arrive as float32
    # and stay float32 through the percentile and zero filters.  Subtracting the
    # cutoffs promotes to float64, because the cutoffs are Python floats, and the
    # published files carry exactly that mix: .X and the "raw" layer float32, the
    # "positives" layer float64.  Casting here would change every stored value.
    X = adata.X
    if sparse.issparse(X):
        X = X.toarray()
    X = np.asarray(X)

    if percentile:
        per_channel = np.percentile(X, p, axis=0)
        keep = (X <= per_channel).all(axis=1)
        adata._inplace_subset_obs(keep)
        X = X[keep, :]

    cutoffs = None
    if subtract_fmo_cutoff or add_fmo_cutoff:
        cutoffs = np.array([marker_cutoff[name] for name in adata.var_names])

    if subtract_fmo_cutoff:
        X = X - cutoffs

    if zero_qc:
        keep = (X > 0).any(axis=1)
        adata._inplace_subset_obs(keep)
        X = X[keep, :]

    if add_fmo_cutoff:
        X = X + cutoffs

    if negatives_to_zeroes:
        X[X < 0] = 0

    if shift:
        X = X - X.min(axis=0)

    if scaling:
        max_per_channel = X.max(axis=0)
        max_per_channel[max_per_channel == 0] = 1
        X = X / max_per_channel

    adata.X = X


def logicle_transform(adata, t, m, w, a, verbose=True):
    """Logicle-transform ``adata.X`` in place, one channel at a time.

    `t`, `m`, `w` and `a` are lists with one entry per channel, in the order of
    ``adata.var_names``.  Working channel by channel keeps only one column in
    memory at a time, which matters on matrices of tens of millions of cells.

    Returns None.  `adata` is modified in place.
    """
    n_channels = adata.n_vars
    for label, values in (("t", t), ("m", m), ("w", w), ("a", a)):
        if len(values) != n_channels:
            raise ValueError(
                f"logicle parameter {label!r} has {len(values)} entries but "
                f"adata has {n_channels} channels"
            )

    if sparse.issparse(adata.X):
        adata.X = adata.X.toarray()

    for index in range(n_channels):
        transformed = pm.tl.normalize_logicle(
            adata[:, index],
            t=t[index],
            m=m[index],
            w=w[index],
            a=a[index],
            inplace=False,
        )
        adata.X[:, index] = transformed.X[:, 0]
        if verbose:
            print(f"channel {index} ({adata.var_names[index]}) done")


# =====================================================================
# Embedding
# =====================================================================


def compute_umap(adata, genes_to_exclude=None, n_pcs=50, n_neighbors=15, random_state=0):
    """Return a copy of `adata` carrying principal components, a graph and a UMAP.

    The embedding is computed on a z-scored copy of the intensities, because
    `scanpy.pp.scale` is applied before the principal component analysis.  The
    returned object keeps the untransformed ``.X`` of the input, so plotting a
    marker still shows the logicle intensity rather than a z-score.  Only
    ``obsm['X_pca']``, ``obsm['X_umap']``, ``obsp['distances']``,
    ``obsp['connectivities']``, ``uns['neighbors']`` and ``uns['umap']`` are
    copied across.

    `genes_to_exclude` names channels to drop before scaling, for example a
    lineage-tracing channel that is not a surface marker.

    `random_state` is passed to the principal component analysis and to UMAP.
    Its default of 0 is the scanpy default, so leaving it alone reproduces the
    published embeddings.  UMAP is only reproducible for a fixed random state
    together with fixed versions of `umap-learn` and `pynndescent`.
    """
    adata_new = adata.copy()

    if genes_to_exclude is not None:
        mask = ~adata.var_names.isin(genes_to_exclude)
        adata_filtered = adata[:, mask].copy()
    else:
        adata_filtered = adata.copy()

    sc.pp.scale(adata_filtered)
    sc.tl.pca(adata_filtered, n_comps=n_pcs, random_state=random_state)
    sc.pp.neighbors(adata_filtered, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata_filtered, random_state=random_state)

    adata_new.obsm["X_pca"] = adata_filtered.obsm["X_pca"]
    adata_new.obsm["X_umap"] = adata_filtered.obsm["X_umap"]
    adata_new.obsp["distances"] = adata_filtered.obsp["distances"]
    adata_new.obsp["connectivities"] = adata_filtered.obsp["connectivities"]
    adata_new.uns["neighbors"] = adata_filtered.uns["neighbors"]
    adata_new.uns["umap"] = adata_filtered.uns["umap"]

    return adata_new


def rotate_umap(adata, angle=-180):
    """Rotate ``adata.obsm['X_umap']`` by `angle` degrees about its centroid.

    Rotation only changes how the embedding is drawn.  It changes no distance,
    no neighbour and no cluster label.  Its purpose is to orient the published
    figures the same way across loops.

    Returns None.  `adata` is modified in place.
    """
    umap_original = adata.obsm["X_umap"].copy()
    center = umap_original.mean(axis=0)
    centred = umap_original - center
    theta = np.deg2rad(angle)
    rotation = np.array(
        [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
    )
    adata.obsm["X_umap"] = centred @ rotation + center


# =====================================================================
# Plotting
# =====================================================================


def plot_marker(adata, markers, threshold=2, n_cols=6, layer=None):
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
            raise ValueError("adata.raw is not set; either set it or pass layer=...")
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
        sb.histplot(np.asarray(X[:, i]).ravel(), bins=bins, ax=axes[i])
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
