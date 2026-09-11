"""Cleaning and logicle transformation of spectral flow cytometry intensities.

`clean_intensities` gathers the optional per-cell and per-channel cleaning steps
that the notebooks switch on and off.  `logicle_transform` applies the logicle
transformation one channel at a time, which is what keeps peak memory down on a
matrix of tens of millions of cells.
"""

import pickle
import json

import numpy as np
import pytometry as pm
from scipy import sparse

__all__ = [
    "load_logicle_parameters",
    "load_marker_cutoffs",
    "clean_intensities",
    "logicle_transform",
]


def load_logicle_parameters(path="logicle_parameters.json"):
    """Return the logicle parameters `(t, m, w, a)` written by the FMO notebook.

    Each of the four is a list with one entry per channel, in the channel order
    of the AnnData object that produced them.  `t` is the top of the display
    range, `m` the number of decades, `w` the number of decades compressed into
    the linear region near zero, and `a` the number of additional negative
    decades.
    """
    with open(path) as handle:
        data = json.load(handle)
    return data["t"], data["m"], data["w"], data["a"]


def load_marker_cutoffs(path="marker_cutoff.pkl"):
    """Return the per-antibody background cutoff dictionary from the FMO notebook.

    A cutoff is the highest intensity reached by the fluorescence-minus-one
    control for that antibody after outlier removal, so an intensity above it is
    treated as real signal.
    """
    with open(path, "rb") as handle:
        return pickle.load(handle)


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
    is set.

    1. `percentile`: drop a cell when any channel is above the `p`-th percentile
       of that channel.  This removes the extreme upper tail of the intensity
       distribution.  Cells are removed, so ``adata.obs`` shrinks with ``adata.X``.
    2. `subtract_fmo_cutoff`: subtract the per-antibody cutoff from every
       intensity, which puts the fluorescence-minus-one background at zero.
    3. `zero_qc`: drop a cell when no channel is above zero.  Run this after the
       subtraction to drop cells that carry no signal in any marker.
    4. `add_fmo_cutoff`: add the per-antibody cutoff back on.  Only meaningful
       together with `subtract_fmo_cutoff`, and then only in combination with a
       step between the two that depends on the shifted values, such as `zero_qc`.
    5. `negatives_to_zeroes`: replace every negative intensity with zero.
    6. `shift`: subtract the per-channel minimum, moving each channel's minimum
       to zero.
    7. `scaling`: divide each channel by its maximum, putting it on a zero to one
       scale.  A channel whose maximum is zero is left alone.

    `marker_cutoff` maps an antibody name to its cutoff, and is only read when
    `subtract_fmo_cutoff` or `add_fmo_cutoff` is set.  Every name in
    ``adata.var_names`` must be a key of it, or a `KeyError` is raised.

    Returns None.  `adata` is modified in place.
    """
    X = adata.X
    if sparse.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

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
    for name in ("t", "m", "w", "a"):
        values = {"t": t, "m": m, "w": w, "a": a}[name]
        if len(values) != n_channels:
            raise ValueError(
                f"logicle parameter {name!r} has {len(values)} entries but "
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
