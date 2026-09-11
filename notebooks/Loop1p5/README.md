# Loop1p5

Plate 2 of experiment LPHO013, eleven experiments in columns 1 to 11.

## Pipeline

| stage | notebook | what it does | reads | writes |
|---|---|---|---|---|
| 1 | `000a_fcs_concat_exp_labelling.ipynb` | read every well, attach its experiment number and culture conditions, concatenate | the `.fcs` wells under `$LABCOMPASS_FCS_ROOT/LPHO013/plate2`<br>`metadata.csv`<br>`live.npz`<br>`counts.npz` | `Loop1p5.h5ad`<br>`Loop1p5_500k.h5ad` |
| 2 | `002a_gating_apply.ipynb` | apply the three gates to every cell | `Loop1p5.h5ad`<br>`gates.npz` | `Loop1p5_Gated.h5ad`<br>`Loop1p5_Gated_500k.h5ad` |
| 3 | `004_transf_singlecell.ipynb` | clean, build the `positives` view, logicle-transform, downsample, embed | `Loop1p5_Gated.h5ad`<br>`marker_cutoff.pkl`<br>`logicle_parameters.json` | `Loop1p5_Raw_Pos.h5ad`<br>`Loop1p5_logicle.h5ad`<br>`Loop1p5_logicle_500k.h5ad`<br>`Loop1p5_UMAP.h5ad` |

Every `.h5ad` file listed under "writes" is produced in this folder and is read
by the stage below it, so the stages run in the order given and only in that
order.

```
000a_fcs_concat_exp_labelling.ipynb
    reads  the .fcs wells under $LABCOMPASS_FCS_ROOT/LPHO013/plate2, metadata.csv, live.npz, counts.npz
    writes Loop1p5.h5ad, Loop1p5_500k.h5ad
002a_gating_apply.ipynb
    reads  Loop1p5.h5ad, gates.npz
    writes Loop1p5_Gated.h5ad, Loop1p5_Gated_500k.h5ad
004_transf_singlecell.ipynb
    reads  Loop1p5_Gated.h5ad, marker_cutoff.pkl, logicle_parameters.json
    writes Loop1p5_Raw_Pos.h5ad, Loop1p5_logicle.h5ad, Loop1p5_logicle_500k.h5ad, Loop1p5_UMAP.h5ad
```

## Experiments in this loop

11 experiments: 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238

The experiment numbers are set in the notebook rather than read from the
`.fcs` files, because the cytometer does not record them.  How a well is
mapped to one of these numbers is described in `000a`.

## Files committed in this folder

- `metadata.csv` — culture conditions, one row per experiment
- `gates.npz` — the three gate polygons, drawn in Loop1 and copied here unchanged
- `live.npz` — the same three gates written against the raw channel names, used to count live cells per well
- `counts.npz` — two gates that select the counting beads, used to count beads per well
- `marker_cutoff.pkl` — one background cutoff per antibody, made in Loop0 by `003` and copied here unchanged
- `logicle_parameters.json` — the logicle parameters `t`, `m`, `w`, `a`, made in Loop0 by `003` and copied here unchanged

Everything else in the table is an intermediate `.h5ad` file.  Those are not
committed, because they run to tens of gigabytes.  Run the notebooks in order to
rebuild them.

## Where the raw data lives

```
/rds/project/rds-SDzz0CATGms/unsorted/for_Juan_Licyel/LPHO013/plate2
```

One unmixed `.fcs` file per well, written by the spectral unmixing software.
The raw data is far too large to commit, so only the path is recorded.  If your
copy sits elsewhere, set the environment variable `LABCOMPASS_FCS_ROOT` before
starting Jupyter rather than editing the notebook.

## A note on file names

The original notebooks wrote this loop's files with the prefix
`BloodPlus_Loop3`, not `Loop1p5`.  That prefix carries another loop's number, because the folder the
notebooks ran in was numbered differently from the folder they were
handed over in.  Reading `BloodPlus_Loop3*.h5ad` as "the files of
that other loop" would be wrong.  The clean notebooks use `Loop1p5`
throughout.  If you have files from the original run, the mapping is
`BloodPlus_Loop3*.h5ad` to `Loop1p5*.h5ad`, and the contents are the same.

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
