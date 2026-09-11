# Loop2p5

Four further experiments of LPHO014.

## Pipeline

| stage | notebook | what it does | reads | writes |
|---|---|---|---|---|
| 1 | `000a_fcs_concat_exp_labelling.ipynb` | read every well, attach its experiment number and culture conditions, concatenate | the `.fcs` wells under `$LABCOMPASS_FCS_ROOT/LPHO014`<br>`metadata_loop2p5.csv`<br>`wells.txt`<br>`live.npz`<br>`counts.npz` | `Loop2p5.h5ad`<br>`Loop2p5_500k.h5ad` |
| 2 | `002a_gating_apply.ipynb` | apply the three gates to every cell | `Loop2p5.h5ad`<br>`gates.npz` | `Loop2p5_Gated.h5ad`<br>`Loop2p5_Gated_500k.h5ad` |
| 3 | `004_transf_singlecell.ipynb` | clean, build the `positives` view, logicle-transform, downsample, embed | `Loop2p5_Gated.h5ad`<br>`marker_cutoff.pkl`<br>`logicle_parameters.json` | `Loop2p5_Raw_Pos.h5ad`<br>`Loop2p5_logicle.h5ad`<br>`Loop2p5_logicle_500k.h5ad`<br>`Loop2p5_UMAP.h5ad` |

Every `.h5ad` file listed under "writes" is produced in this folder and is read
by the stage below it, so the stages run in the order given and only in that
order.

```
000a_fcs_concat_exp_labelling.ipynb
    reads  the .fcs wells under $LABCOMPASS_FCS_ROOT/LPHO014, metadata_loop2p5.csv, wells.txt, live.npz, counts.npz
    writes Loop2p5.h5ad, Loop2p5_500k.h5ad
002a_gating_apply.ipynb
    reads  Loop2p5.h5ad, gates.npz
    writes Loop2p5_Gated.h5ad, Loop2p5_Gated_500k.h5ad
004_transf_singlecell.ipynb
    reads  Loop2p5_Gated.h5ad, marker_cutoff.pkl, logicle_parameters.json
    writes Loop2p5_Raw_Pos.h5ad, Loop2p5_logicle.h5ad, Loop2p5_logicle_500k.h5ad, Loop2p5_UMAP.h5ad
```

## Experiments in this loop

4 experiments: 248, 249, 250, 251

The experiment numbers are set in the notebook rather than read from the
`.fcs` files, because the cytometer does not record them.  How a well is
mapped to one of these numbers is described in `000a`.

## Files committed in this folder

- `metadata_loop2p5.csv` — culture conditions, one row per experiment
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
/rds/project/rds-SDzz0CATGms/unsorted/for_Juan_Licyel/LPHO014
```

One unmixed `.fcs` file per well, written by the spectral unmixing software.
The raw data is far too large to commit, so only the path is recorded.  If your
copy sits elsewhere, set the environment variable `LABCOMPASS_FCS_ROOT` before
starting Jupyter rather than editing the notebook.

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
