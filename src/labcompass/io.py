"""Reading spectral flow cytometry wells and attaching per-well metadata.

Each well of a plate is one ``.fcs`` file.  The functions here read a well,
attach the metadata that the cytometer wrote into the file, attach the
experimental conditions from a lookup table, and finally concatenate every well
into one AnnData object.

How a well is mapped to an experiment number differs between loops, because the
plate layouts differ, so that mapping stays in each loop's own notebook.  What
is shared is everything that happens once the experiment number is known.
"""

import glob
import os

import anndata as ann
import numpy as np
import pandas as pd
import pytometry as pm

__all__ = [
    "CHANNELS_TO_DROP",
    "LIVE_DEAD_CHANNEL",
    "find_fcs_files",
    "load_metadata_table",
    "attach_cytometer_metadata",
    "attach_condition_metadata",
    "concat_wells",
    "add_live_dead_obs",
    "drop_channels",
    "rename_channels_to_antibody",
]


# Channels that carry no antibody signal and are dropped before analysis.  The
# seven "[AF color N]" channels are the unmixed autofluorescence components that
# the spectral cytometer reports.  The live/dead channel is moved to ``.obs``
# by `add_live_dead_obs` before being dropped here.  CD184 was dropped because
# its staining failed in this panel.
CHANNELS_TO_DROP = [
    "CD184 : PerCP-eFluor710 - Area",
    "Live_Dead : 7-AAD - Area",
    "[AF color 1] - Area",
    "[AF color 2] - Area",
    "[AF color 3] - Area",
    "[AF color 4] - Area",
    "[AF color 5] - Area",
    "[AF color 6] - Area",
    "[AF color 7] - Area",
]

LIVE_DEAD_CHANNEL = "Live_Dead : 7-AAD - Area"


def find_fcs_files(data_path, pattern="*AutoSpectral.fcs", subdir="filtered"):
    """Return the sorted list of unmixed ``.fcs`` files under `data_path`.

    `subdir` restricts the search to directories with that name, which is where
    the unmixing software writes its output.  Pass None to search every
    directory below `data_path`.

    The list is sorted, so the order of the wells does not depend on the order
    in which the filesystem happens to return them.
    """
    parts = [data_path, "**"]
    if subdir is not None:
        parts.append(subdir)
    parts.append(pattern)
    files = glob.glob(os.path.join(*parts), recursive=True)
    files.sort()
    return files


def load_metadata_table(path, sep="\t", index_col="experiment_number"):
    """Read the per-experiment condition table and return it indexed by string.

    The index is forced to `str`.  Without that, pandas infers an integer index
    when every experiment number happens to be numeric, and a later lookup with
    a string key silently finds nothing.  Forcing the type makes the lookup
    behave the same whichever experiment numbers a loop happens to use.
    """
    metadata = pd.read_csv(path, sep=sep)
    metadata = metadata.set_index(index_col)
    metadata.index = metadata.index.astype(str)
    return metadata


def attach_cytometer_metadata(adata):
    """Copy the acquisition metadata from ``adata.uns['meta']`` into ``adata.obs``.

    Four columns are added: ``date``, ``cytometer``, ``cytometer_serial_no`` and
    ``well_id``.  Every cell of the well gets the same value, because these
    describe the acquisition rather than the cell.
    """
    meta = adata.uns["meta"]
    adata.obs["date"] = meta["date"]
    adata.obs["cytometer"] = meta["cyt"]
    adata.obs["cytometer_serial_no"] = meta["cytsn"]
    adata.obs["well_id"] = meta["smno"]


def attach_condition_metadata(adata, metadata, experiment_number):
    """Copy the row of `metadata` for `experiment_number` into ``adata.obs``.

    One column is added per column of `metadata`, holding the culture conditions
    of that experiment such as cytokine concentrations.  `experiment_number` is
    converted to `str` to match the index type that `load_metadata_table` sets.

    A `KeyError` is raised when the experiment number is not in the table, which
    is deliberate: silently skipping the row would give cells with no recorded
    condition.
    """
    key = str(experiment_number)
    if key not in metadata.index:
        raise KeyError(
            f"experiment_number {key!r} is not in the metadata table; "
            f"the table holds {len(metadata)} experiments"
        )
    adata.obs[metadata.columns.tolist()] = metadata.loc[key].tolist()


def concat_wells(adatas):
    """Concatenate a list of per-well AnnData objects into one object.

    Only the channels present in every well are kept, by ``join="inner"``.  Cell
    names are made unique across wells by appending the well's position in the
    list, by ``index_unique="-"``, which is what lets later steps index cells by
    name without collisions.
    """
    return ann.concat(
        adatas,
        join="inner",
        merge="same",
        uns_merge="unique",
        label="batch",
        index_unique="-",
    )


def add_live_dead_obs(adata, channel=LIVE_DEAD_CHANNEL, obs_name="7-AAD"):
    """Copy the live/dead channel into ``adata.obs[obs_name]``.

    The gating strategy draws its live gate against this value, and gating reads
    from ``.obs``, so the channel has to be copied across before it is dropped
    from ``.X`` by `drop_channels`.
    """
    values = adata[:, channel].X
    if hasattr(values, "toarray"):
        values = values.toarray()
    adata.obs[obs_name] = np.asarray(values).ravel()


def drop_channels(adata, channels_to_drop=None):
    """Return `adata` without the named channels.

    Defaults to `CHANNELS_TO_DROP`.  Names that are not present are ignored, so
    a panel that never had a given autofluorescence channel does not fail.
    """
    if channels_to_drop is None:
        channels_to_drop = CHANNELS_TO_DROP
    keep = ~np.isin(adata.var["marker"].values, list(channels_to_drop))
    return adata[:, keep]


def rename_channels_to_antibody(adata):
    """Rename the channels from ``"CD34 : APC - Area"`` to ``"CD34"``, in place.

    A channel name from the unmixing software has the form
    ``"<antibody> : <fluorochrome> - <parameter>"``.  This splits it, writes the
    antibody name into ``adata.var["antibody"]`` and into ``adata.var_names``,
    and returns the list of fluorochromes in channel order.

    A `ValueError` is raised when two channels reduce to the same antibody name,
    because duplicate variable names break later lookups by name.
    """
    markers = list(adata.var["marker"].values)
    antibodies = [marker.split(" :")[0] for marker in markers]
    fluorochromes = [marker.split(": ")[1].split(" -")[0] for marker in markers]

    duplicates = sorted({a for a in antibodies if antibodies.count(a) > 1})
    if duplicates:
        raise ValueError(f"channels reduce to duplicate antibody names: {duplicates}")

    rename = dict(zip(markers, antibodies))

    adata.var = adata.var.copy()
    adata.var["antibody"] = (
        pd.Categorical(adata.var["marker"]).map(rename, na_action=None).astype("category")
    )
    adata.var_names = antibodies
    return fluorochromes


def read_well(path):
    """Read one ``.fcs`` file and return it with string cell names and marker names.

    ``adata.var_names`` is set to the marker column so that a gate saved against
    a channel name can be applied before ``pytometry.pp.split_signal`` moves the
    scatter channels into ``.obs``.
    """
    adata = pm.io.read_fcs(path)
    adata.obs.index = adata.obs.index.astype(str)
    adata.var_names = adata.var["marker"].astype(str)
    return adata
