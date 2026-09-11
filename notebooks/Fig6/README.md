# Fig6

Experiment LPHO015, the twenty experiments behind Figure 6.  The pipeline is the same as for the numbered loops; only the input experiment and the downstream figure panels differ.

## Pipeline

| stage | notebook | what it does | reads | writes |
|---|---|---|---|---|
| 1 | `000a_fcs_concat_exp_labelling.ipynb` | read every well, attach its experiment number and culture conditions, concatenate | the `.fcs` wells under `$LABCOMPASS_FCS_ROOT/LPHO015`<br>`metadata_results.csv`<br>`wells.txt`<br>`live.npz`<br>`counts.npz` | `Fig6.h5ad`<br>`Fig6_500k.h5ad` |
| 2 | `002a_gating_apply.ipynb` | apply the three gates to every cell | `Fig6.h5ad`<br>`gates.npz` | `Fig6_Gated.h5ad`<br>`Fig6_Gated_500k.h5ad` |
| 3 | `004_transf_singlecell.ipynb` | clean, build the `positives` view, logicle-transform, downsample, embed, cluster | `Fig6_Gated.h5ad`<br>`marker_cutoff.pkl`<br>`logicle_parameters.json` | `Fig6_Raw_Pos.h5ad`<br>`Fig6_logicle.h5ad`<br>`Fig6_logicle_500k.h5ad`<br>`Fig6_UMAP.h5ad`<br>`Fig6_Leiden.h5ad` |

Every `.h5ad` file listed under "writes" is produced in this folder and is read
by the stage below it, so the stages run in the order given and only in that
order.

```
000a_fcs_concat_exp_labelling.ipynb
    reads  the .fcs wells under $LABCOMPASS_FCS_ROOT/LPHO015, metadata_results.csv, wells.txt, live.npz, counts.npz
    writes Fig6.h5ad, Fig6_500k.h5ad
002a_gating_apply.ipynb
    reads  Fig6.h5ad, gates.npz
    writes Fig6_Gated.h5ad, Fig6_Gated_500k.h5ad
004_transf_singlecell.ipynb
    reads  Fig6_Gated.h5ad, marker_cutoff.pkl, logicle_parameters.json
    writes Fig6_Raw_Pos.h5ad, Fig6_logicle.h5ad, Fig6_logicle_500k.h5ad, Fig6_UMAP.h5ad, Fig6_Leiden.h5ad
```

## Experiments in this loop

20 experiments: 262, 263, 264, 265, 266, 267, 270, 271, 273, 276, 275, 261, 279, 277, 280, 269, 268, 272, 278, 274

The experiment numbers are set in the notebook rather than read from the
`.fcs` files, because the cytometer does not record them.  How a well is
mapped to one of these numbers is described in `000a`.

## Files committed in this folder

- `metadata_results.csv` — culture conditions, one row per experiment
- `wells.txt` — the wells to read, in the order that decides which experiment each well belongs to
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
/rds/project/rds-SDzz0CATGms/unsorted/for_Juan_Licyel/LPHO015
```

One unmixed `.fcs` file per well, written by the spectral unmixing software.
The raw data is far too large to commit, so only the path is recorded.  If your
copy sits elsewhere, set the environment variable `LABCOMPASS_FCS_ROOT` before
starting Jupyter rather than editing the notebook.

## A note on file names

The original notebooks wrote this loop's files with the prefix
`results`, not `Fig6`.  That prefix says nothing about which loop the file belongs to.  The clean notebooks use `Fig6`
throughout.  If you have files from the original run, the mapping is
`results*.h5ad` to `Fig6*.h5ad`, and the contents are the same.

## Parameters

| parameter | value | set in |
|---|---|---|
| cells kept per experiment after downsampling | about 500,000 in total | `000a`, `002a`, `004` |
| upper percentile filter | 99.9997 | `004` |
| principal components | 19 | `004` |
| neighbours for the graph and UMAP | 15 | `004` |
| random state for the principal components and UMAP | 0 | `004` |
| Leiden resolution | 0.9 | `004` |
| random seed for Leiden | 51 | `004` |
| UMAP rotation, in degrees | 45 | `004` |
| random seed for downsampling | 0 | every stage that downsamples |
