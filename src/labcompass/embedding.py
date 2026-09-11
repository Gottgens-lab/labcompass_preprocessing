"""Principal components, neighbourhood graph and UMAP for gated cytometry data."""

import numpy as np
import scanpy as sc

__all__ = ["compute_umap", "rotate_umap"]


def compute_umap(
    adata,
    genes_to_exclude=None,
    n_pcs=50,
    n_neighbors=15,
    random_state=0,
):
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
    published embeddings.  Note that UMAP is only reproducible for a fixed
    random state together with a fixed version of ``umap-learn`` and of
    ``pynndescent``.
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
