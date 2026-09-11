"""Balanced downsampling of an AnnData object by an observation column."""

import anndata as ann
import numpy as np

__all__ = ["downsample", "target_size_for"]


def target_size_for(adata, desired_cell_number, obs="experiment_number"):
    """Return the per-group cap that gives roughly `desired_cell_number` cells.

    The cap is the integer part of `desired_cell_number` divided by the number
    of distinct values of `adata.obs[obs]`.  The total after downsampling is
    smaller than `desired_cell_number` whenever a group holds fewer cells than
    the cap, because such a group is kept whole rather than padded.
    """
    n_groups = len(adata.obs[obs].unique())
    return int(desired_cell_number / n_groups)


def downsample(
    adata: ann.AnnData,
    obs: str = "experiment_replicate",
    target_size: int = None,
    seed: int = 0,
) -> ann.AnnData:
    """Cap every group in ``adata.obs[obs]`` at `target_size` cells.

    A group with more than `target_size` cells is sampled without replacement.
    A group with fewer cells is kept whole, so the groups are only balanced when
    every group is at least `target_size` cells large.  When `target_size` is
    None the size of the smallest group is used, which does balance the groups.

    Returns a new AnnData object.  The input is not modified.

    Reproducibility: this function calls ``np.random.seed(seed)``, which reseeds
    the global NumPy generator as a side effect.  Two calls with the same seed,
    the same group sizes and the same NumPy version select the same cells.  Any
    code that draws random numbers after this call is affected by the reseeding.
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

    keep_positions = np.concatenate(keep_positions)
    return adata[keep_positions].copy()
