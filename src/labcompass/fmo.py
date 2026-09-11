"""Background cutoffs and logicle parameters from fluorescence-minus-one controls.

A fluorescence-minus-one control is a sample stained with every antibody of the
panel except one.  Whatever signal appears in the missing antibody's channel is
therefore background, so the top of that signal marks where real staining
begins.  `compute_marker_cutoffs` measures that top for every antibody and turns
it into the parameters that the logicle transformation needs.
"""

import json
import math
import pickle

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb
from scipy.stats import gaussian_kde

__all__ = ["compute_marker_cutoffs", "save_marker_cutoffs", "save_logicle_parameters"]


def compute_marker_cutoffs(
    adata,
    antibody_list,
    fmo_obs="fmo",
    fully_stained_label="fully_stained",
    fmo_percentile=99.5,
    stained_percentile=99.75,
    plot=True,
):
    """Return the background cutoff and the logicle parameters for each antibody.

    For antibody ``A`` the function takes the cells whose ``adata.obs[fmo_obs]``
    equals ``"A"``, which are the cells of the control that lacks ``A``, reads
    channel ``A`` for those cells, drops everything above the `fmo_percentile`-th
    percentile as outlier noise, and takes the maximum of what is left.  That
    maximum is the cutoff.

    Four lists are returned alongside, in the order of `antibody_list`, holding
    the logicle parameters:

    - ``t``: the maximum of channel ``A`` among the fully stained cells, after
      dropping everything above the `stained_percentile`-th percentile.  This is
      the top of the display range.
    - ``m``: ``log10(t)``, the number of decades that the display spans.
    - ``w``: ``log10(cutoff)``, the number of decades compressed into the linear
      region around zero.
    - ``a``: ``-log10(cutoff)``, the number of additional negative decades shown.

    Returns ``(marker_cutoff, t, m, w, a)`` where `marker_cutoff` maps an
    antibody name to its cutoff.

    An antibody with no control cells, or whose surviving maximum is not
    positive, raises a `ValueError` rather than being skipped.  Skipping would
    leave the four lists shorter than `antibody_list` and silently misalign
    every later channel.
    """
    marker_cutoff = {}
    t, m, w, a = [], [], [], []

    stained_mask = (adata.obs[fmo_obs] == fully_stained_label).values
    if not stained_mask.any():
        raise ValueError(
            f"no cells with {fmo_obs} == {fully_stained_label!r}; "
            "the fully stained sample has to be concatenated in first"
        )
    stained = adata[stained_mask, :].copy()

    if plot:
        fig, axes = plt.subplots(
            nrows=len(antibody_list), figsize=(6, len(antibody_list) * 3)
        )
        axes = np.atleast_1d(axes)
    else:
        fig, axes = None, [None] * len(antibody_list)

    for i, marker in enumerate(antibody_list):
        control_mask = (adata.obs[fmo_obs] == marker).values
        if not control_mask.any():
            raise ValueError(f"no control cells for antibody {marker!r}")

        channel_mask = (adata.var["antibody"] == marker).values
        control = _dense_column(adata[control_mask, :][:, channel_mask])

        cutoff = np.percentile(control, fmo_percentile)
        control_clean = control[control <= cutoff]
        if control_clean.size == 0:
            raise ValueError(f"no control cells left for {marker!r} after outlier removal")

        top_of_background = float(control_clean.max())
        if top_of_background <= 0:
            raise ValueError(
                f"the control maximum for {marker!r} is {top_of_background}, "
                "which has no logarithm; check the channel assignment"
            )

        stained_values = _dense_column(stained[:, channel_mask])
        stained_cutoff = np.percentile(stained_values, stained_percentile)
        stained_clean = stained_values[stained_values <= stained_cutoff]
        top_of_scale = float(stained_clean.max())
        if top_of_scale <= 0:
            raise ValueError(
                f"the stained maximum for {marker!r} is {top_of_scale}, "
                "which has no logarithm"
            )

        marker_cutoff[marker] = top_of_background
        t.append(top_of_scale)
        m.append(math.log10(top_of_scale))
        w.append(math.log10(top_of_background))
        a.append(-math.log10(top_of_background))

        if plot:
            _plot_one(axes[i], control_clean, stained_clean, marker)

    if plot:
        plt.tight_layout()

    return marker_cutoff, t, m, w, a


def _dense_column(view):
    """Return the single-channel intensities of `view` as a flat float array."""
    X = view.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.asarray(X).ravel()


def _plot_one(ax, control_clean, stained_clean, marker):
    """Draw the control and the fully stained distribution of one antibody."""
    xmin = -3 * control_clean.max()
    xmax = 6 * control_clean.max()
    bins = np.linspace(xmin, xmax, 41)

    sb.histplot(control_clean, bins=bins, ax=ax, stat="density", color="red")
    sb.histplot(stained_clean, bins=bins, ax=ax, stat="density", color="black")

    kde = gaussian_kde(control_clean)
    x_eval = np.linspace(xmin, xmax, 500)
    ax.plot(x_eval, kde(x_eval), color="red", label="control density")

    ax.set_title(f"FMO = {marker}, background reaches {control_clean.max():.1f}")
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel("intensity")
    ax.set_ylabel("density")


def save_marker_cutoffs(marker_cutoff, path="marker_cutoff.pkl"):
    """Write the antibody-to-cutoff dictionary to `path` as a pickle."""
    with open(path, "wb") as handle:
        pickle.dump(marker_cutoff, handle)


def save_logicle_parameters(t, m, w, a, path="logicle_parameters.json"):
    """Write the four logicle parameter lists to `path` as JSON."""
    with open(path, "w") as handle:
        json.dump({"t": t, "m": m, "w": w, "a": a}, handle, indent=2)
