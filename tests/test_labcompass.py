"""Tests for `notebooks/labcompass.py`.

They run on a small synthetic AnnData object, so they finish in seconds and need
no access to the raw data.  Run them from the repository root with:

    pytest tests -v
"""

import anndata as ann
import numpy as np
import pandas as pd
import pytest

from labcompass import (
    apply_saved_gates,
    clean_intensities,
    compute_umap,
    downsample,
    get_axis_data,
    plot_ab_histograms,
    plot_leiden_clusters,
    plot_marker,
    plot_umap_by_group,
    read_gate_file,
    rotate_umap,
)

N_CELLS = 2000
N_CHANNELS = 6
MARKERS = [f"M{i}" for i in range(N_CHANNELS)]


@pytest.fixture
def adata():
    """A small object shaped like a concatenated plate: six markers, five experiments."""
    rng = np.random.default_rng(0)
    var = pd.DataFrame(
        {
            "marker": [f"M{i} : FL{i} - Area" for i in range(N_CHANNELS)],
            "antibody": pd.Categorical(MARKERS),
        },
        index=MARKERS,
    )
    obs = pd.DataFrame(
        {
            "experiment_number": rng.integers(1, 6, N_CELLS).astype(str),
            "FSC-A :: FSC - Area": rng.uniform(0, 1e6, N_CELLS),
            "SSC-A :: SSC - Area": rng.uniform(0, 1e6, N_CELLS),
        },
        index=[f"c{i}" for i in range(N_CELLS)],
    )
    return ann.AnnData(X=np.abs(rng.normal(500, 200, (N_CELLS, N_CHANNELS))), obs=obs, var=var)


def write_gates(path, polygons, x_cols, y_cols):
    np.savez(
        path,
        polygons=np.array(polygons, dtype=object),
        x_cols=np.array(x_cols, dtype=object),
        y_cols=np.array(y_cols, dtype=object),
    )


# ------------------------------------------------------------------ gating
def test_get_axis_data_finds_obs_and_var(adata):
    assert get_axis_data(adata, "FSC-A :: FSC - Area").shape == (N_CELLS,)
    assert get_axis_data(adata, "M0").shape == (N_CELLS,)
    with pytest.raises(KeyError):
        get_axis_data(adata, "not a channel")


def test_apply_saved_gates_combines_with_and(adata, tmp_path):
    """mode="all" keeps a cell only when every polygon contains it."""
    path = str(tmp_path / "gates.npz")
    everything = np.array([[0, 0], [1e6, 0], [1e6, 1e6], [0, 1e6]], float)
    left_half = np.array([[0, 0], [5e5, 0], [5e5, 1e6], [0, 1e6]], float)
    write_gates(path, [everything, left_half],
                ["FSC-A :: FSC - Area"] * 2, ["SSC-A :: SSC - Area"] * 2)

    gated, masks = apply_saved_gates(adata, path, mode="all")
    assert masks[0].sum() == N_CELLS
    assert gated.n_obs == masks[1].sum() < N_CELLS
    assert (gated.obs["FSC-A :: FSC - Area"] <= 5e5).all()

    gated_any, _ = apply_saved_gates(adata, path, mode="any")
    assert gated_any.n_obs == N_CELLS


def test_gate_file_round_trip_with_equal_vertex_counts(tmp_path):
    """numpy stores equal-length polygons as one 3-D array; reading must still work."""
    path = str(tmp_path / "gates.npz")
    square = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
    write_gates(path, [square, square], ["a", "a"], ["b", "b"])
    polygons, x_cols, y_cols = read_gate_file(path)
    assert len(polygons) == len(x_cols) == len(y_cols) == 2
    assert polygons[0].shape == (4, 2)


def test_gate_file_round_trip_with_unequal_vertex_counts(tmp_path):
    path = str(tmp_path / "gates.npz")
    square = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
    triangle = np.array([[0, 0], [1, 0], [0, 1]], float)
    write_gates(path, [square, triangle], ["a", "a"], ["b", "b"])
    polygons, _, _ = read_gate_file(path)
    assert [p.shape for p in polygons] == [(4, 2), (3, 2)]


def test_empty_gate_file_is_rejected(adata, tmp_path):
    path = str(tmp_path / "gates.npz")
    write_gates(path, [], [], [])
    with pytest.raises(ValueError, match="no polygons"):
        apply_saved_gates(adata, path)


def test_gate_on_a_marker_channel(adata, tmp_path):
    """A gate may name a marker channel, not only an .obs column.

    This is what `live.npz` and `counts.npz` need, and it is the one behaviour
    that the two versions of `apply_saved_gates` in the original notebooks
    disagreed on.
    """
    path = str(tmp_path / "gates.npz")
    box = np.array([[0, 0], [600, 0], [600, 600], [0, 600]], float)
    write_gates(path, [box], ["M0"], ["M1"])
    gated, _ = apply_saved_gates(adata, path)
    assert 0 < gated.n_obs < N_CELLS


def test_gate_with_an_unknown_axis_fails_loudly(adata, tmp_path):
    path = str(tmp_path / "gates.npz")
    box = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
    write_gates(path, [box], ["no such axis"], ["M1"])
    with pytest.raises(KeyError):
        apply_saved_gates(adata, path)


# ---------------------------------------------------------------- sampling
def test_downsample_caps_each_group(adata):
    small = downsample(adata, obs="experiment_number", target_size=200, seed=0)
    assert (small.obs["experiment_number"].value_counts() <= 200).all()


def test_downsample_is_reproducible(adata):
    a = downsample(adata, obs="experiment_number", target_size=50, seed=0)
    b = downsample(adata, obs="experiment_number", target_size=50, seed=0)
    assert list(a.obs_names) == list(b.obs_names)


def test_downsample_keeps_small_groups_whole(adata):
    """A group below the cap is kept whole rather than padded."""
    adata.obs["experiment_number"] = ["1"] * 10 + ["2"] * (N_CELLS - 10)
    small = downsample(adata, obs="experiment_number", target_size=100, seed=0)
    counts = small.obs["experiment_number"].value_counts()
    assert counts["1"] == 10
    assert counts["2"] == 100


def test_downsample_selects_by_position_not_by_label(adata):
    """Duplicate cell names must not pull in extra rows."""
    adata.obs_names = ["dup"] * N_CELLS
    small = downsample(adata, obs="experiment_number", target_size=50, seed=0)
    assert small.n_obs == small.obs["experiment_number"].value_counts().sum()
    assert (small.obs["experiment_number"].value_counts() <= 50).all()


# --------------------------------------------------------------- transform
def test_clean_intensities_keeps_obs_and_X_aligned(adata):
    cutoffs = {name: 400.0 for name in adata.var_names}
    clean_intensities(adata, cutoffs, percentile=True, p=99.0, zero_qc=True)
    assert adata.X.shape[0] == adata.obs.shape[0] == adata.n_obs


def test_clean_intensities_percentile_drops_the_upper_tail(adata):
    before = adata.n_obs
    clean_intensities(adata, {}, percentile=True, p=99.0)
    assert adata.n_obs < before


def test_clean_intensities_scaling_gives_zero_to_one(adata):
    cutoffs = {name: 400.0 for name in adata.var_names}
    clean_intensities(
        adata, cutoffs, subtract_fmo_cutoff=True, negatives_to_zeroes=True, scaling=True
    )
    assert adata.X.min() >= 0.0
    assert adata.X.max() == pytest.approx(1.0)


def test_clean_intensities_keeps_float32_as_float32(adata):
    """The filters must not promote the dtype.

    The published files carry .X and the "raw" layer as float32 and the
    "positives" layer as float64.  That mix is not arbitrary: subtracting the
    cutoffs promotes to float64 because the cutoffs are Python floats, while the
    percentile and zero filters only select rows.  Forcing a dtype anywhere in
    this function changes every stored value.
    """
    adata.X = adata.X.astype(np.float32)
    clean_intensities(adata, {}, percentile=True, p=99.0, zero_qc=True)
    assert adata.X.dtype == np.float32

    cutoffs = {name: 400.0 for name in adata.var_names}
    clean_intensities(adata, cutoffs, subtract_fmo_cutoff=True)
    assert adata.X.dtype == np.float64


def test_clean_intensities_rejects_an_unknown_antibody(adata):
    with pytest.raises(KeyError):
        clean_intensities(adata, {"not a marker": 1.0}, subtract_fmo_cutoff=True)


def test_clean_intensities_subtract_then_add_is_a_round_trip(adata):
    """add_fmo_cutoff undoes subtract_fmo_cutoff when nothing between them acts."""
    before = adata.X.copy()
    cutoffs = {name: 400.0 for name in adata.var_names}
    clean_intensities(adata, cutoffs, subtract_fmo_cutoff=True, add_fmo_cutoff=True)
    assert np.allclose(before, adata.X)


# --------------------------------------------------------------- embedding
def test_compute_umap_leaves_X_untouched(adata):
    small = downsample(adata, obs="experiment_number", target_size=60, seed=0)
    embedded = compute_umap(small, n_pcs=3, n_neighbors=10, random_state=0)
    assert embedded.obsm["X_umap"].shape == (small.n_obs, 2)
    assert np.allclose(embedded.X, small.X)


def test_rotate_umap_is_a_rotation(adata):
    """A full turn returns every point to where it started."""
    adata.obsm["X_umap"] = np.random.default_rng(0).normal(5, 3, (N_CELLS, 2))
    before = adata.obsm["X_umap"].copy()
    rotate_umap(adata, angle=90)
    assert not np.allclose(before, adata.obsm["X_umap"])
    rotate_umap(adata, angle=270)
    assert np.allclose(before, adata.obsm["X_umap"], atol=1e-8)


def test_rotate_umap_preserves_distances(adata):
    """Rotation must not move any cell relative to any other."""
    adata = adata[:200].copy()
    adata.obsm["X_umap"] = np.random.default_rng(0).normal(5, 3, (200, 2))
    from scipy.spatial.distance import pdist

    before = pdist(adata.obsm["X_umap"])
    rotate_umap(adata, angle=37)
    assert np.allclose(before, pdist(adata.obsm["X_umap"]), atol=1e-9)


# ---------------------------------------------------------------- plotting
def test_plot_helpers_return_a_figure(adata):
    import matplotlib

    matplotlib.use("Agg")
    adata.obsm["X_umap"] = np.random.default_rng(0).normal(size=(N_CELLS, 2))
    adata.raw = ann.AnnData(X=adata.X, var=adata.var.copy(), obs=adata.obs.copy())
    adata.obs["leiden"] = pd.Categorical(
        np.random.default_rng(0).integers(0, 4, N_CELLS).astype(str)
    )

    assert plot_marker(adata, markers=["M0", "M1"], threshold=400) is not None
    assert plot_leiden_clusters(adata) is not None
    assert plot_ab_histograms(adata, bins=10) is not None
    assert plot_umap_by_group(adata, obs="experiment_number") is not None


def test_plot_marker_rejects_an_unknown_marker(adata):
    import matplotlib

    matplotlib.use("Agg")
    adata.obsm["X_umap"] = np.zeros((N_CELLS, 2))
    adata.raw = ann.AnnData(X=adata.X, var=adata.var.copy(), obs=adata.obs.copy())
    with pytest.raises(KeyError):
        plot_marker(adata, markers=["not a marker"])
