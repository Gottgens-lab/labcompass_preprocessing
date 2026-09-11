# What the original notebooks did, and what changed

This repository is a cleaned rewrite of the notebooks at
`/rds/project/rds-SDzz0CATGms/users/jar82/collaborations/theis/Notebooks`.
This file records what was checked, what was found, and what was changed, so
that a reader can tell a rewrite from a correction.

Every claim below was checked against the original notebook files. The cell
numbers are zero-based indices into the `cells` list of the `.ipynb` file, which
is what a script sees; Jupyter shows the same cells without numbering markdown.

## Summary

No error was found that changes a published number. Five defects would have
produced a wrong result if the notebooks had been rerun from a clean kernel, and
they are listed under "Defects that would fire on a rerun". The rest are
robustness and clarity problems.

The largest single problem is not a bug but a habit: the notebooks were run out
of order. 11 of the 32 notebooks store execution counts that do not increase
down the file, and 8 store the same execution count in more than one cell. A
notebook in that state does not say what it did, because the state a cell saw
depends on the order the person clicked, not on the order the cells sit in.

## Where the original code lived

The `Notebooks/` folder is a copy made for handover. The notebooks were run from
a set of sibling folders one level higher, which is why their relative data path
`../../../../../unsorted/for_Juan_Licyel/` resolves correctly from those folders
and not from `Notebooks/<Loop>/`. The mapping is:

| folder here | `Notebooks/` folder | original working folder | original file prefix |
|---|---|---|---|
| `Loop0` | `Loop0` | `theis/BloodPlus` | `BloodPlus_AS` |
| `Loop1` | `Loop1` | `theis/BloodPlus_Loop2` | `BloodPlus_Loop2` |
| `Loop1p5` | `Loop1p5` | `theis/BloodPlus_Loop3` | `BloodPlus_Loop3` |
| `Loop2` | `Loop2` | `theis/Loop2` | `Loop2` |
| `Loop2p5` | `Loop2p5` | `theis/Loop2p5` | `Loop2p5` |
| `Loop3` | `Loop3` | `theis/Loop3` | `Loop3` |
| `Loop4` | `Loop4` | `theis/Loop4` | `Loop4` |
| `Loop4p5` | `Loop4p5` | `theis/Loop4p5` | `Loop4p5` |
| `Fig6` | `Fig6` | `theis/Results_Fig6` | `results` |

Two of those prefixes name a different loop than the folder they sit in. The
files that `Loop1` wrote are called `BloodPlus_Loop2*.h5ad`, and the files that
`Loop1p5` wrote are called `BloodPlus_Loop3*.h5ad`. The prefix came from the
original working folder name, which was numbered differently. The clean
notebooks use the loop's own name as the prefix everywhere.

## Defects that would fire on a rerun

### 1. Loop0 stage 04: the 5 million cell subset is drawn from the 500 thousand cell subset

In `Loop0/004_transf_singlecell.ipynb`, cell 32 downsamples to about 500,000
cells and rebinds `adata` to the result. Cell 35 then downsamples again, to
about 5,000,000 cells, and writes the result as `BloodPlus_Logicle_pm_5M.h5ad`.
Read top to bottom, cell 35 downsamples an object that already holds only
500,000 cells, so the file named `_5M` would hold 500,000 cells.

The stored execution counts show this did not happen in the recorded run: cell
32 has count 62 and cell 35 has count 16, so cell 35 ran first, while `adata`
still held all 52,085,765 cells. The file on disk is therefore correct and the
notebook is wrong.

Fixed by keeping the full object under one name and drawing every subset from
it, so the order of the cells and the order of execution agree.

### 2. Loop0 stage 04: one cluster's cells would silently become missing values

Cell 64 maps 25 Leiden clusters to cell type names in `annotation_dict`, then
fixes the display order with

```python
adata.obs['cell_type'] = pd.Categorical(adata.obs['cell_type'], categories=new_order, ordered=True)
```

`new_order` holds 25 entries, but it omits `LMPP-GMPMono`, which
`annotation_dict` assigns to cluster 22, and it includes an empty string `''`
that names no cluster. Any cell in cluster 22 would become a missing value with
no warning, because `pd.Categorical` drops a value that is not in `categories`.

The recorded run produced 21 clusters, numbered 0 to 20, so cluster 22 did not
exist and the defect did not fire. It would fire on any rerun that produced 23
or more clusters.

Fixed by adding `LMPP-GMPMono` to the order, removing the empty string, and
adding a check that raises if any name in the mapping is missing from the order
or if any cluster has no name.

### 3. Loop0 stage 04: a cell reads a file that nothing writes

Cell 56 is `adata = sc.read("BloodPlus_pm_UMAP_500k_.h5ad")`. No cell in any
notebook of the project writes that name; the nearest is cell 54, which writes
`BloodPlus_pm_UMAP_5M.h5ad`. Note the trailing underscore in the name that is
read. The cell has no stored execution count, so it never ran in the recorded
session.

Fixed by removing the reread. The clean stage 04 keeps one object in memory
through the embedding and the clustering.

### 4. Loop0 stage 04: an indentation error

Cell 62 is ` plot_leiden_clusters(adata)`, with one leading space. That is an
`IndentationError` in Python. The cell has a stored execution count of 4, which
means it did run, so the leading space was added after that run.

Fixed by removing the space.

### 5. Loop0 stage 01b: two functions are called but never defined

In `Loop0/001b_fcs_concat_FMOs.ipynb`, cell 7 calls `is_int(x)` and
`is_float(x)` inside a loop over the `meta_`-prefixed keys of `adata.uns["meta"]`.
Neither function is defined anywhere in the project, which was checked by
searching every notebook. The calls would raise `NameError`.

They did not fire because the control files have no `meta_` keys, so the list
`positions` is empty and the loop body never runs. The cell completed with
execution count 6.

Fixed by removing the block. The control wells carry no experimental metadata,
so nothing was lost.

## Defects in the shared helper functions

The same block of helper functions was pasted into the first code cell of nearly
every notebook. Four distinct versions of that block exist, identified by
checksum. All 32 copies are byte-identical to one of the four.

### 6. Two versions of `apply_saved_gates` that are not interchangeable

The version in the `001a` and `002a` notebooks reads a gate axis only from
`.obs`:

```python
pts = np.c_[adata.obs[x].values, adata.obs[y].values]
```

The version in the `000a` notebooks reads it from `.obs` or from `.var_names`,
through a helper called `get_axis_data`. Only the second version can apply the
gates in `live.npz` and `counts.npz`, whose axis names are channel names such as
`'Live_Dead : 7-AAD - Area'` and `'CD235a : AF700 - Area'`. Calling the first
version on those files raises `KeyError`.

The notebooks happen to pair each version with gate files it can handle, so the
defect never fired. It is still a trap: moving a cell between notebooks changes
which version is in scope.

Fixed by keeping one version, the one that reads from either place. It is a
strict superset in behaviour, so no gating result changes.

### 7. `add_fmo_CO` on its own does nothing

In `transform`, the cutoff is subtracted whenever either flag is set:

```python
if subtract_fmo_CO or add_fmo_CO:
    cutoffs = np.array([marker_cutoff[gene] for gene in adata.var_names])
    X = X - cutoffs
...
if add_fmo_CO:
    X += cutoffs
```

Calling with `add_fmo_CO=True` and `subtract_fmo_CO=False` therefore subtracts
the cutoff and adds it straight back, which changes nothing. Every call in every
notebook passes `add_fmo_CO=False`, so no published number depends on this.

Changed in the rewrite: `clean_intensities` subtracts only when
`subtract_fmo_cutoff` is set, and adds only when `add_fmo_cutoff` is set. That
is a behaviour change for a combination that the notebooks never use.

### 8. Duplicate definitions

In the `000a` block, `get_axis_data` and `downsample` are each defined twice,
one after the other, with identical bodies. Harmless, and removed.

### 9. Dead code with a missing import

Four functions were never called from any notebook: `sweep_umap`,
`umap_experiments`, `plot_meld_density` and `plot_umap_kde_per_condition`.
`plot_meld_density` calls `meld.MELD(...)`, and `meld` is neither imported in any
notebook nor installed in the analysis environment, so the function could not
run as written.

Removed from the package. If they are wanted back, they are in the original
notebooks under the same names.

### 10. `downsample` reseeds the global random generator

`downsample` calls `np.random.seed(seed)`, which reseeds NumPy's global
generator as a side effect. Any random draw made later in the same session is
affected. The behaviour was kept, because changing it would select different
cells and therefore change the published subsets, but it is now documented in
the docstring.

The original also selected cells by their `.obs` index labels. The rewrite
selects by position instead, which avoids a wrong result if two cells ever share
a name. That the two give the same cells for the same seed was checked directly:
`np.random.choice(labels, k, replace=False)` and
`labels[np.random.choice(len(labels), k, replace=False)]` return the same
elements after the same seeding.

## Things that are fragile rather than wrong

### 11. The experiment lookup depended on an inferred pandas dtype

Loop0 stage 00 tests `experiment_list[i]['experiment_number'] in metadata.index`,
comparing a string read from a CSV against the index of a table read with
`pd.read_csv`. If pandas had inferred an integer index, every test would have
failed and every well would have been skipped without an error. It worked
because the index happened to be read as text.

Fixed: `load_metadata_table` forces the index to `str`.

### 12. A well that is skipped leaves no trace

In every stage 00, a well whose experiment number is not in the condition table
is skipped by an `if` with no `else`. Nothing records that it happened, so a
mismatch between the well list and the table shows up only as a smaller cell
count.

Fixed: skipped wells are collected and printed with the reason.

### 13. The plate-group mapping can run off the end of the list

Loop2, Loop2p5, Loop3, Loop4, Loop4p5 and Fig6 assign experiment numbers by
counting groups: each new combination of plate name and plate row takes the next
entry of a hard-coded list. If the wells produce more groups than the list has
entries, `experiment_list[expno]` raises `IndexError` at an arbitrary point,
after some wells have already been labelled.

Fixed: the clean notebook checks the bound before indexing and raises a message
that names the well and both counts.

### 14. The well list was rebuilt three times in Loop2

`Loop2/000a_fcs_concat_exp_labelling.ipynb` builds `files` by pattern in cell 5,
filters five wells out of it in cells 6 and 7, then overwrites it with a
hard-coded list in cell 8, then overwrites that with a shorter hard-coded list in
cell 9. Only the list in cell 9 has any effect; cells 5 to 8 are dead.

The order of the well list is not cosmetic. It is what decides which experiment
number each well is given, under the group-counting rule above.

Fixed: the effective list is stored in `wells.txt` next to the notebook, with
paths relative to the data root, and the notebook reads it and checks that every
file exists. The dead cells are gone.

### 15. Commented-out code left in place

Loop0 stage 02 carries the earlier version of `experiments_to_remove` in a
triple-quoted string above the version that runs, and the two lists differ.
Loop0 stage 04 cell 25 holds a logicle implementation inside a string, ending
with the comment `#Needs to be fixed`, while cell 27 holds the version that ran.
Loop0 stage 03 cells 17 and 18 hold usage examples in strings, one of which
names a file, `my_lists.json`, that does not exist.

All removed. The version that ran is the version that is kept.

### 16. A typo that happened to be harmless

Loop0 stage 03 cell 13 reads:

```python
percentile = 99.5 #might work with 99.
5
```

The `5` of a second `99.5` ended up on its own line, where Python evaluates it
and throws it away. The percentile actually used is 99.5. Removed.

### 17. `seed` was set but not used for the embedding

Stage 04 sets `seed = 51` and passes it to `sc.tl.leiden`, but `compute_umap`
calls `sc.tl.pca` and `sc.tl.umap` without a random state, so both used the
scanpy default of 0. The published embeddings were therefore made with random
state 0, and `seed` only ever affected the clustering.

In the rewrite `compute_umap` takes `random_state` explicitly, and the clean
stage 04 sets `EMBEDDING_RANDOM_STATE = 0`, which is what reproduces the
published embeddings. The Leiden seed is separate and stays 51, under the name
`LEIDEN_SEED`. Nothing moves; the two seeds are just no longer confused with
each other.

### 17b. Only two loops clustered

`sc.tl.leiden` and `rotate_umap` appear only in `Loop0` (resolution 0.8) and
`Fig6` (resolution 0.9). The other seven loops stop after the UMAP and a panel
grid showing where each experiment sits in it.

The clean notebooks follow that: the seven loops that never clustered do not
cluster now. Adding a clustering step to a loop that never ran one would put
un-reviewed results into the repository.

### 18. `counter.pkl` and the partial-save machinery

Stage 04 wrote a channel counter to `counter.pkl` so that a logicle run killed
part way through could restart. The flag that enables it, `temporalsave`, is
`False` in every notebook, so the counter is always 0 and the partial save never
happens. Removed, along with the file.

## Fig6 depends on a dataset that is not in this repository

`Fig6/004_transf_singlecell.ipynb` cell 57 is

```python
adata = sc.read("ds_HID01_logicle_100k_UMAP_ann.h5ad")
```

That file is produced by no notebook in the project and is not in the folder.
From that cell onward the notebook works on that dataset rather than on Fig6's
own output, including the cell type annotation and the frame-by-frame UMAP
figures.

The clean `Fig6` folder covers stages 00, 02 and 04 up to the Leiden clustering,
which is the part that runs from data this repository has. The steps that need
`ds_HID01_logicle_100k_UMAP_ann.h5ad` are not included, because a notebook that
cannot run is not publishable. Add that file and its provenance, and those steps
can be restored.

## What was checked and found correct

- The 20-antibody panel, the channels dropped, and the antibody renaming are
  identical across all nine loops.
- `marker_cutoff.pkl` and `logicle_parameters.json` are byte-identical in all
  nine folders, so every loop uses the same intensity calibration. Checked by
  checksum.
- `gates.npz` has two distinct versions: Loop0 has its own, and the other eight
  folders share the one drawn in Loop1. Checked by checksum.
- The cleaning and embedding parameters are the same in every loop: percentile
  99.9997, 19 principal components, 15 neighbours, Leiden resolution 0.8, UMAP
  rotation 45 degrees. Fig6 uses Leiden resolution 0.9.
- Every well named in every `wells.txt` exists under the data root. Checked by
  listing the filesystem: 100 of 100 for Fig6, 52 of 52 for Loop2, 42 of 42 for
  Loop3, 21 of 21 for Loop2p5, 6 of 6 for Loop4 and Loop4p5.
- The three gate polygons in `gates.npz` are drawn on forward scatter against
  side scatter, forward scatter area against forward scatter height, and forward
  scatter against 7-AAD, in that order. Read directly from the file.

## How the rewrite was checked

- Every code cell of all 32 clean notebooks parses as Python.
- No clean notebook stores an output or an execution count.
- Every name that a clean notebook imports from `labcompass` exists in the
  package. Checked by importing the package and looking each name up.
- The library has 20 tests in `tests/`, all passing. They cover the gate file
  round trip for both storage shapes that numpy produces, gating on an `.obs`
  column and on a marker channel, the `and` and `or` combination rules,
  downsampling for cap, reproducibility and small groups, the alignment of
  `.obs` with `.X` after cleaning, the zero-to-one scaling, the rotation of the
  embedding, and the error paths.

What was **not** checked: no clean notebook has been run end to end on the real
data. Doing so needs the analysis environment, which is not installed here, and
the logicle step alone takes several hours per loop. The clean notebooks are
therefore written and reviewed, not executed.
