"""Preprocessing for the LabCompass spectral flow cytometry screen.

The pipeline turns one ``.fcs`` file per well into one AnnData object per
experimental loop, carrying logicle-transformed marker intensities, a UMAP
embedding and Leiden clusters.  It runs in five numbered stages, which the
notebooks under ``notebooks/<Loop>/`` follow one file each.  See the repository
README for the stage-by-stage description and each loop's README for the exact
files that loop reads and writes.

Import the pieces you need, for example::

    from labcompass import apply_saved_gates, downsample, clean_intensities
"""

__version__ = "1.0.0"

from .embedding import compute_umap, rotate_umap
from .fmo import (
    compute_marker_cutoffs,
    save_logicle_parameters,
    save_marker_cutoffs,
)
from .gating import (
    DrawGate,
    apply_saved_gates,
    get_axis_data,
    init_gate_file,
    read_gate_file,
)
from .io import (
    CHANNELS_TO_DROP,
    LIVE_DEAD_CHANNEL,
    add_live_dead_obs,
    attach_condition_metadata,
    attach_cytometer_metadata,
    concat_wells,
    drop_channels,
    find_fcs_files,
    load_metadata_table,
    read_well,
    rename_channels_to_antibody,
)
from .plotting import (
    plot_ab_histograms,
    plot_leiden_clusters,
    plot_marker,
    plot_umap_by_group,
)
from .sampling import downsample, target_size_for
from .transform import (
    clean_intensities,
    load_logicle_parameters,
    load_marker_cutoffs,
    logicle_transform,
)

__all__ = [
    "__version__",
    # io
    "CHANNELS_TO_DROP",
    "LIVE_DEAD_CHANNEL",
    "find_fcs_files",
    "read_well",
    "load_metadata_table",
    "attach_cytometer_metadata",
    "attach_condition_metadata",
    "concat_wells",
    "add_live_dead_obs",
    "drop_channels",
    "rename_channels_to_antibody",
    # gating
    "DrawGate",
    "get_axis_data",
    "init_gate_file",
    "read_gate_file",
    "apply_saved_gates",
    # fmo
    "compute_marker_cutoffs",
    "save_marker_cutoffs",
    "save_logicle_parameters",
    # transform
    "load_logicle_parameters",
    "load_marker_cutoffs",
    "clean_intensities",
    "logicle_transform",
    # sampling
    "downsample",
    "target_size_for",
    # embedding
    "compute_umap",
    "rotate_umap",
    # plotting
    "plot_marker",
    "plot_leiden_clusters",
    "plot_ab_histograms",
    "plot_umap_by_group",
]
