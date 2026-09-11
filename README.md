# labcompass_preprocessing

Preprocessing for the LabCompass spectral flow cytometry screen: from one `.fcs`
file per well to one annotated `AnnData` object per experimental loop, carrying
logicle-transformed marker intensities, a UMAP embedding and Leiden clusters.

The screen measures a 20-antibody haematopoietic panel on cord blood cells grown
under many combinations of cytokines and small molecules.  Each combination is
one *experiment*; each round of the screen is one *loop*.

## What is in here

```
notebooks/labcompass.py   the few helpers that are too long to repeat in every notebook
notebooks/<Loop>/         one folder per loop: its notebooks, its lookup tables, its README
env/                      the pinned environment that produced the published files
docs/AUDIT.md             what was found in the original notebooks and what was changed
tests/                    tests for notebooks/labcompass.py
```

There is nothing to install.  `labcompass.py` sits beside the loop folders and a
notebook reaches it with `sys.path.append("..")`.  It holds only the pieces that
are long and shared: the interactive gating class, gate application,
downsampling, the cleaning steps, the logicle loop, the embedding and the four
plotting functions.  Everything else a stage does is written out in the notebook
cell that does it, so you can read a stage without opening a second file.

Each loop folder has its own `README.md` giving that loop's file-dependency
graph: which notebook reads which file and which notebook writes it.  Read that
first when you want to rerun a loop.

## The loops

| loop | what it covers | experiments | raw data |
|---|---|---|---|
| `Loop0` | every plate acquired up to that point; draws the gates and builds the calibration | 122 after removals | all of `for_Juan_Licyel/` |
| `Loop1` | LPHO013 plate 1 | 216 to 227 | `LPHO013/plate1` |
| `Loop1p5` | LPHO013 plate 2 | 228 to 238 | `LPHO013/plate2` |
| `Loop2` | LPHO014, first block | 239 to 247 | `LPHO014` |
| `Loop2p5` | LPHO014, second block | 248 to 251 | `LPHO014` |
| `Loop3` | LPHO014, third block | 252 to 258 | `LPHO014` |
| `Loop4` | LPHO014, fourth block | 259 | `LPHO014` |
| `Loop4p5` | LPHO014, fifth block | 260 | `LPHO014` |
| `Fig6` | LPHO015, the experiments behind Figure 6 | 261 to 280 | `LPHO015` |

`Loop0` is the reference loop.  It is the only one that builds the
fluorescence-minus-one calibration, and with `Loop1` one of only two that draw
gates by hand.  Only `Loop0` and `Fig6` go on to cluster; the other seven stop
after the UMAP.  Every
other loop carries a copy of `marker_cutoff.pkl` and `logicle_parameters.json`
and uses them unchanged, so the intensity scale is the same across loops.

## The pipeline

Five stages.  Not every loop runs every stage; the loop's README says which.

| stage | notebook | what it does |
|---|---|---|
| 00 | `000a_fcs_concat_exp_labelling.ipynb` | read each well, attach its experiment number and culture conditions, concatenate into one object |
| 01 | `001a_gating_maker.ipynb` | draw the three gates by hand (`Loop0` and `Loop1` only) |
| 01b | `001b_fcs_concat_FMOs.ipynb` | read the fluorescence-minus-one control wells (`Loop0` only) |
| 02 | `002a_gating_apply.ipynb` | apply the gates to every cell |
| 02b | `002b_gating_apply_FMOs.ipynb` | apply the same gates to the controls (`Loop0` only) |
| 03 | `003_FMO_maker.ipynb` | measure the background per antibody, derive the logicle parameters (`Loop0` only) |
| 04 | `004_transf_singlecell.ipynb` | clean, transform, downsample, embed, and cluster in `Loop0` and `Fig6` |

### What each stage means

**Gating.** Three polygons, applied in order, each on the cells the previous one
kept.  Forward scatter area against side scatter area separates cells from
debris.  Forward scatter area against forward scatter height removes doublets,
because two cells passing the laser together give a wide pulse.  Forward scatter
area against the 7-AAD channel removes dead cells, because 7-AAD only enters a
cell whose membrane is broken.

**The fluorescence-minus-one calibration.** A fluorescence-minus-one control is
a sample stained with every antibody of the panel except one.  Whatever signal
appears in the missing antibody's channel is background rather than staining, so
the top of that signal marks where real staining begins.  Stage 03 measures that
top for each of the 20 antibodies and writes it to `marker_cutoff.pkl`.

**The logicle transformation.** Linear near zero and logarithmic further out.
Spectral unmixing produces negative values, which a plain logarithm cannot show,
so logicle is what keeps the negative and near-zero events visible.  Its four
parameters `t`, `m`, `w` and `a` are derived from the same controls and stored
in `logicle_parameters.json`.

**The `positives` layer.** A second view of the same cells, in which each
channel has had its background cutoff subtracted, negative values set to zero,
and the rest scaled to run from zero to one.  A value above zero here means the
marker is really present, which the logicle values cannot say on their own,
because the logicle scale has no fixed zero point.  The marker panels in stage
04 read this layer through `adata.raw`.

## Getting started

```bash
conda env create -f env/environment.yml
conda activate labcompass
jupyter lab
```

Then open the loop you want under `notebooks/` and run its notebooks in
numerical order.

### Where the raw data lives

```
/rds/project/rds-SDzz0CATGms/unsorted/for_Juan_Licyel
```

One unmixed `.fcs` file per well.  The raw data is far too large to commit, so
the repository records only the path.  If your copy sits elsewhere, set the
environment variable before starting Jupyter:

```bash
export LABCOMPASS_FCS_ROOT=/your/path/to/for_Juan_Licyel
```

Every notebook reads that variable and falls back to the path above.

### What is committed and what is not

Committed: the notebooks, `labcompass.py`, the per-experiment condition
tables (`metadata*.csv`), the well lists (`wells.txt`), the gate polygons
(`gates.npz`, `live.npz`, `counts.npz`) and the calibration files
(`marker_cutoff.pkl`, `logicle_parameters.json`).  These total under 400 KB and
are what makes the published results reproducible without the raw data.

Not committed: every `.h5ad` file and every `.fcs` file.  The full objects run to
tens of gigabytes.  `.gitignore` excludes them.

## Reproducibility

This was checked against the original outputs, not only asserted.  Three loops
were rebuilt from their raw `.fcs` wells and every file compared with the one the
original notebooks produced.  All 14 files agree exactly on every intensity:
`.X`, both layers, `.raw`, and all 46 `.obs` columns.  See `docs/AUDIT.md` for
the file-by-file table.

The one thing that does not carry across machines is the UMAP embedding.  It is
reproducible on a given machine with the pinned environment, and two runs there
agree to the last bit, but it differs from the published embedding by a median of
0.122 UMAP units.  UMAP optimises its layout with a parallel stochastic gradient
descent whose update order depends on the thread count, so the coordinates depend
on the machine as well as on the seed.  Everything upstream of the embedding is
portable.  If you need the published coordinates, read them from the published
`.h5ad` rather than recomputing them.

Every notebook runs top to bottom.  No cell reads a file that a later cell
writes, so `Restart Kernel and Run All Cells` rebuilds that notebook's outputs
from its inputs.  The one exception is `001a_gating_maker.ipynb`, which opens an
interactive figure and waits for a person; you do not need to run it, because
the gates it produces are committed.

Three things fix the numbers:

- Three separate seeds, each set at the top of the notebook that uses it.
  Downsampling uses 0.  The principal components and UMAP use random state 0,
  which is the scanpy default and is what the published embeddings were made
  with.  Leiden uses 51.  In the original notebooks a single variable called
  `seed` was set to 51 and reached Leiden but never reached UMAP, which made it
  look as though the embeddings used 51 when they did not.
- The gate polygons and the calibration files are committed, so the gating and
  the intensity scale do not have to be redrawn or remeasured.
- The package versions are pinned in `env/requirements.txt`.  UMAP and Leiden
  are only reproducible for a fixed seed together with fixed versions of
  `umap-learn`, `pynndescent`, `numba` and `leidenalg`.  Changing any of those
  four moves the embedding and renumbers the clusters, which in turn invalidates
  any mapping from cluster number to cell type.

The long step is the logicle transformation in stage 04.  On `Loop4p5`, at
642,482 cells, it takes about 10 seconds per channel for 20 channels.  On
`Loop0`, at 52 million cells, it takes about 10 minutes per channel.

## Relation to the original notebooks

This repository is a cleaned rewrite.  `docs/AUDIT.md` lists what was found in
the original notebooks, what was fixed, and what was deliberately left alone.
It also gives the mapping from the original working directories and file names
to the ones used here.
