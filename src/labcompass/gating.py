"""Polygon gating for spectral flow cytometry data.

Two things live here:

`DrawGate` opens an interactive matplotlib figure in which the user drags the
vertices of a polygon over a two-channel scatter plot, then appends the polygon
to a ``.npz`` gate file.  It needs an interactive matplotlib backend
(``%matplotlib widget`` plus ``ipympl``) and therefore only works inside Jupyter.

`apply_saved_gates` replays a saved gate file on any AnnData object without a
display, so the batch part of the pipeline is fully non-interactive.

A gate file is a ``.npz`` archive with three parallel arrays.  Entry ``i`` of
``polygons`` is an ``(n_vertices, 2)`` array of vertex coordinates, and entries
``i`` of ``x_cols`` and ``y_cols`` name the two axes that the polygon was drawn
against.  An axis name is looked up first in ``adata.obs`` and then in
``adata.var_names``, so a gate may be drawn against a scatter channel such as
``"FSC-A :: FSC - Area"``, a derived observation such as ``"7-AAD"``, or a
marker channel such as ``"CD235a : AF700 - Area"``.
"""

from collections import Counter, defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pytometry as pm
import scanpy as sc
from matplotlib.patches import Polygon
from matplotlib.path import Path
from scipy.spatial import ConvexHull, Delaunay

__all__ = [
    "get_axis_data",
    "DrawGate",
    "apply_saved_gates",
    "init_gate_file",
    "read_gate_file",
]


def get_axis_data(adata, key):
    """Return a one-dimensional float array for `key`, from `.obs` or from `.X`.

    `key` is looked up in ``adata.obs.columns`` first and in ``adata.var_names``
    second.  A `KeyError` is raised when it is in neither, which is what makes a
    typo in a gate file fail loudly instead of producing an empty gate.
    """
    if key in adata.obs.columns:
        return adata.obs[key].values

    if key in adata.var_names:
        vals = adata[:, key].X
        if hasattr(vals, "toarray"):
            vals = vals.toarray()
        return np.asarray(vals).ravel()

    raise KeyError(f"{key!r} not found in adata.obs or adata.var_names")


def init_gate_file(filename):
    """Create an empty gate file, discarding any gate file already at `filename`.

    Call this once before drawing the first gate of a gating strategy.
    `DrawGate.save_gate` appends, so skipping this step silently keeps the gates
    of the previous run.
    """
    np.savez(
        filename,
        polygons=np.array([], dtype=object),
        x_cols=np.array([], dtype=object),
        y_cols=np.array([], dtype=object),
    )


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
    which is the sequential gating strategy that the notebooks use.  With
    ``mode="any"`` a cell is kept when it falls inside at least one polygon.

    Returns a tuple of the gated copy of `adata` and the list of per-gate boolean
    masks, each mask having one entry per cell of the input `adata`.
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

    # ------------------------------------------------------------------
    # mouse interaction
    # ------------------------------------------------------------------

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

    # ------------------------------------------------------------------
    # results
    # ------------------------------------------------------------------

    def get_masks(self):
        """Return one boolean mask per polygon, over the cells of `self.adata`."""
        pts = np.c_[
            get_axis_data(self.adata, self.x), get_axis_data(self.adata, self.y)
        ]
        return [Path(poly).contains_points(pts) for poly in self.verts]

    def filter_cells(self, mode="any"):
        """Return the cells that this gate keeps, as a new AnnData object."""
        masks = self.get_masks()
        if mode == "any":
            final = np.logical_or.reduce(masks)
        else:
            final = np.logical_and.reduce(masks)
        return self.adata[final].copy()

    def save_gate(self, filename):
        """Append this gate's polygon and its two axis names to `filename`.

        The gate file is created when it does not exist, so calling
        `init_gate_file` first is what keeps a re-run from stacking new gates on
        top of the gates of the previous run.
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

    # ------------------------------------------------------------------
    # automatic starting polygons
    # ------------------------------------------------------------------

    @staticmethod
    def _bounding_rect(pts):
        """Return the four vertices of the axis-aligned bounding box of `pts`."""
        pts = np.asarray(pts)
        xmin, ymin = pts.min(axis=0)
        xmax, ymax = pts.max(axis=0)
        return np.array(
            [[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]]
        )

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
            if nxt == start:
                break
            if nxt in visited:
                break
            ring.append(nxt)
            visited.add(nxt)
            prev, current = current, nxt

        if len(ring) < 3:
            return pts[ConvexHull(pts).vertices]

        return pts[np.array(ring, dtype=int)]


# The notebooks were written against a lowercase class name.  Keep it working.
drawgate = DrawGate
