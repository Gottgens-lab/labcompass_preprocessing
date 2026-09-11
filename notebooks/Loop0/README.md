# Loop0

The first and largest loop.  It pools every plate acquired up to that point, draws the gating strategy by hand, builds the fluorescence-minus-one calibration that every later loop reuses, and drops the experiments that were not part of the screen.

## Pipeline

| stage | notebook | what it does | reads | writes |
|---|---|---|---|---|
| 1 | `000a_fcs_concat_exp_labelling.ipynb` | read every well, attach its experiment number and culture conditions, concatenate | the `.fcs` wells under `$LABCOMPASS_FCS_ROOT/`<br>`experiment_table.csv`<br>`metadata.csv` | `Loop0.h5ad`<br>`Loop0_500k.h5ad` |
| 2 | `001a_gating_maker.ipynb` | draw the three gates by hand (needs a person, see below) | `Loop0_500k.h5ad` | `gates.npz` |
| 3 | `001b_fcs_concat_FMOs.ipynb` | read the fluorescence-minus-one control wells | the `.fcs` controls under `$LABCOMPASS_FCS_ROOT/LPFMO02` | `Loop0_FMOs.h5ad` |
| 4 | `002a_gating_apply.ipynb` | apply the three gates to every cell, then drop the experiments that are not part of the screen | `Loop0.h5ad`<br>`gates.npz` | `Loop0_Gated.h5ad`<br>`Loop0_Gated_500k.h5ad`<br>`Loop0_Gated_ExpRem.h5ad`<br>`Loop0_Gated_ExpRem_500k.h5ad` |
| 5 | `002b_gating_apply_FMOs.ipynb` | apply the same three gates to the control cells | `Loop0_FMOs.h5ad`<br>`gates.npz` | `Loop0_FMOs_Gated.h5ad` |
| 6 | `003_FMO_maker.ipynb` | measure the background per antibody and derive the logicle parameters | `Loop0_FMOs_Gated.h5ad`<br>`Loop0_Gated_ExpRem.h5ad` | `marker_cutoff.pkl`<br>`logicle_parameters.json` |
| 7 | `004_transf_singlecell.ipynb` | clean, build the `positives` view, logicle-transform, downsample, embed, cluster, name the clusters | `Loop0_Gated_ExpRem.h5ad`<br>`marker_cutoff.pkl`<br>`logicle_parameters.json` | `Loop0_Raw_Pos.h5ad`<br>`Loop0_logicle.h5ad`<br>`Loop0_logicle_500k.h5ad`<br>`Loop0_UMAP.h5ad`<br>`Loop0_Leiden.h5ad`<br>`Loop0_celltype.h5ad` |

Every `.h5ad` file listed under "writes" is produced in this folder and is read
by the stage below it, so the stages run in the order given and only in that
order.

```
000a_fcs_concat_exp_labelling.ipynb
    reads  the .fcs wells under $LABCOMPASS_FCS_ROOT/, experiment_table.csv, metadata.csv
    writes Loop0.h5ad, Loop0_500k.h5ad
001a_gating_maker.ipynb
    reads  Loop0_500k.h5ad
    writes gates.npz
001b_fcs_concat_FMOs.ipynb
    reads  the .fcs controls under $LABCOMPASS_FCS_ROOT/LPFMO02
    writes Loop0_FMOs.h5ad
002a_gating_apply.ipynb
    reads  Loop0.h5ad, gates.npz
    writes Loop0_Gated.h5ad, Loop0_Gated_500k.h5ad, Loop0_Gated_ExpRem.h5ad, Loop0_Gated_ExpRem_500k.h5ad
002b_gating_apply_FMOs.ipynb
    reads  Loop0_FMOs.h5ad, gates.npz
    writes Loop0_FMOs_Gated.h5ad
003_FMO_maker.ipynb
    reads  Loop0_FMOs_Gated.h5ad, Loop0_Gated_ExpRem.h5ad
    writes marker_cutoff.pkl, logicle_parameters.json
004_transf_singlecell.ipynb
    reads  Loop0_Gated_ExpRem.h5ad, marker_cutoff.pkl, logicle_parameters.json
    writes Loop0_Raw_Pos.h5ad, Loop0_logicle.h5ad, Loop0_logicle_500k.h5ad, Loop0_UMAP.h5ad, Loop0_Leiden.h5ad, Loop0_celltype.h5ad
```

## Files committed in this folder

- `experiment_table.csv` — maps an acquisition date and start time to an experiment number
- `metadata.csv` — culture conditions, one row per experiment
- `gates.npz` — the three gate polygons, drawn in this loop by `001a`
- `marker_cutoff.pkl` — one background cutoff per antibody, made here by `003`
- `logicle_parameters.json` — the logicle parameters `t`, `m`, `w`, `a`, made here by `003`

Everything else in the table is an intermediate `.h5ad` file.  Those are not
committed, because they run to tens of gigabytes.  Run the notebooks in order to
rebuild them.

## Where the raw data lives

```
/rds/project/rds-SDzz0CATGms/unsorted/for_Juan_Licyel/
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
`BloodPlus_AS`, not `Loop0`.  For Loop0 that prefix names a
different loop, which was a source of confusion, so the clean notebooks use
`Loop0` throughout.  If you have files from the original run, the mapping is
`BloodPlus_AS*.h5ad` to `Loop0*.h5ad`.  The contents are the same.

## Parameters

| parameter | value | set in |
|---|---|---|
| cells kept per experiment after downsampling | about 500,000 in total | `000a`, `002a`, `004` |
| upper percentile filter | 99.9997 | `004` |
| principal components | 19 | `004` |
| neighbours for the graph and UMAP | 15 | `004` |
| random state for the principal components and UMAP | 0 | `004` |
| Leiden resolution | 0.8 | `004` |
| random seed for Leiden | 51 | `004` |
| UMAP rotation, in degrees | 45 | `004` |
| random seed for downsampling | 0 | every stage that downsamples |
