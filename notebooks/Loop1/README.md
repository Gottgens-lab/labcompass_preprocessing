# Loop1

Plate 1 of experiment LPHO013, twelve experiments in columns 1 to 12.

## Pipeline

| stage | notebook | what it does | reads | writes |
|---|---|---|---|---|
| 1 | `000a_fcs_concat_exp_labelling.ipynb` | read every well, attach its experiment number and culture conditions, concatenate | the `.fcs` wells under `$LABCOMPASS_FCS_ROOT/LPHO013/plate1/filtered`<br>`metadata.csv`<br>`live.npz`<br>`counts.npz` | `Loop1.h5ad`<br>`Loop1_500k.h5ad` |
| 2 | `001a_gating_maker.ipynb` | draw the three gates by hand (needs a person, see below) | `Loop1_500k.h5ad` | `gates.npz` |
| 3 | `002a_gating_apply.ipynb` | apply the three gates to every cell | `Loop1.h5ad`<br>`gates.npz` | `Loop1_Gated.h5ad`<br>`Loop1_Gated_500k.h5ad` |
| 4 | `004_transf_singlecell.ipynb` | clean, build the `positives` view, logicle-transform, downsample, embed | `Loop1_Gated.h5ad`<br>`marker_cutoff.pkl`<br>`logicle_parameters.json` | `Loop1_Raw_Pos.h5ad`<br>`Loop1_logicle.h5ad`<br>`Loop1_logicle_500k.h5ad`<br>`Loop1_UMAP.h5ad` |

Every `.h5ad` file listed under "writes" is produced in this folder and is read
by the stage below it, so the stages run in the order given and only in that
order.

```
000a_fcs_concat_exp_labelling.ipynb
    reads  the .fcs wells under $LABCOMPASS_FCS_ROOT/LPHO013/plate1/filtered, metadata.csv, live.npz, counts.npz
    writes Loop1.h5ad, Loop1_500k.h5ad
001a_gating_maker.ipynb
    reads  Loop1_500k.h5ad
    writes gates.npz
002a_gating_apply.ipynb
    reads  Loop1.h5ad, gates.npz
    writes Loop1_Gated.h5ad, Loop1_Gated_500k.h5ad
004_transf_singlecell.ipynb
    reads  Loop1_Gated.h5ad, marker_cutoff.pkl, logicle_parameters.json
    writes Loop1_Raw_Pos.h5ad, Loop1_logicle.h5ad, Loop1_logicle_500k.h5ad, Loop1_UMAP.h5ad
```

## Experiments in this loop

12 experiments: 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227

The experiment numbers are set in the notebook rather than read from the
`.fcs` files, because the cytometer does not record them.  How a well is
mapped to one of these numbers is described in `000a`.

## Files committed in this folder

- `metadata.csv` — culture conditions, one row per experiment
- `gates.npz` — the three gate polygons, drawn in this loop by `001a`
- `live.npz` — the same three gates written against the raw channel names, used to count live cells per well
- `counts.npz` — two gates that select the counting beads, used to count beads per well
- `marker_cutoff.pkl` — one background cutoff per antibody, made in Loop0 by `003` and copied here unchanged
- `logicle_parameters.json` — the logicle parameters `t`, `m`, `w`, `a`, made in Loop0 by `003` and copied here unchanged

Everything else in the table is an intermediate `.h5ad` file.  Those are not
committed, because they run to tens of gigabytes.  Run the notebooks in order to
rebuild them.

## Where the raw data lives

```
/rds/project/rds-SDzz0CATGms/unsorted/for_Juan_Licyel/LPHO013/plate1/filtered
```

One unmixed `.fcs` file per well, written by the spectral unmixing software.
The raw data is far too large to commit, so only the path is recorded.  If your
copy sits elsewhere, set the environment variable `LABCOMPASS_FCS_ROOT` before
starting Jupyter rather than editing the notebook.

## The one notebook that needs a person

`001a_gating_maker.ipynb` opens an interactive figure and waits for you to drag
polygon vertices.  It cannot run unattended.

You do not need to run it.  `gates.npz` is committed, so `002a` reproduces the
published gating on its own.  Run `001a` only when you want to redraw the gates,
and be aware that redrawing them changes every downstream file.

## A note on file names

The original notebooks wrote this loop's files with the prefix
`BloodPlus_Loop2`, not `Loop1`.  That prefix carries another loop's number, because the folder the
notebooks ran in was numbered differently from the folder they were
handed over in.  Reading `BloodPlus_Loop2*.h5ad` as "the files of
that other loop" would be wrong.  The clean notebooks use `Loop1`
throughout.  If you have files from the original run, the mapping is
`BloodPlus_Loop2*.h5ad` to `Loop1*.h5ad`, and the contents are the same.

## Parameters

| parameter | value | set in |
|---|---|---|
| cells kept per experiment after downsampling | about 500,000 in total | `000a`, `002a`, `004` |
| upper percentile filter | 99.9997 | `004` |
| principal components | 19 | `004` |
| neighbours for the graph and UMAP | 15 | `004` |
| random state for the principal components and UMAP | 0 | `004` |
| Leiden clustering | not run in this loop | |
| UMAP rotation | not applied in this loop | |
| random seed for downsampling | 0 | every stage that downsamples |
