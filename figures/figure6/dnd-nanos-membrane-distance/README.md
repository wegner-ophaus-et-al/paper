# TMD-signal distance from membrane

Analysis of confocal images of zebrafish embryos (6 hpf) to test whether a TMD (transmembrane domain)-tagged construct changes the localization of Dnd1 mRNA/protein "spots" relative to the plasma membrane. For each cell, the pipeline segments the cell and Dnd1 spots, computes each spot's distance from the membrane and the membrane-vs-cytoplasm distribution of Dnd1 signal, and compares these across four injection conditions (`wt`, `tmd-n_stuffer`, `tmd-n_d`, `n_d`).

## Input data

The pipeline expects a root folder of per-sample subfolders, one per imaged cell, laid out as:

```
<root>/<uid>__<condition_string>-<cell_number>/original_raw/<name>.lsm
<root>/<uid>__<condition_string>-<cell_number>/masks/cell.tif      # written/updated by the pipeline
<root>/<uid>__<condition_string>-<cell_number>/masks/spots.tif     # written by the pipeline
```

- `uid` is an 8-hex-character identifier used to blind samples during manual curation.
- `condition_string` is an underscore-separated pair of injection component codes (e.g. `A709_F125`), mapped in `analyze.py` to one of `wt`, `tmd-n_stuffer`, `tmd-n_d`, `n_d`.
- Each `.lsm` file is a 3-channel Zeiss confocal stack: channel 0 = Vasa, channel 1 = Dnd1, channel 2 = membrane.

The raw microscopy data is not included in this repository and is not yet public.

## How to run

No environment/dependency file is present in this repo (no `requirements.txt`/`environment.yml`). The code imports `pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`, `scikit-image`, `tifffile`, `Pillow`, `micro_sam`, and `torch` (with optional `psutil` for memory checks on Apple Silicon); set up an environment with these yourself.

1. Edit the `root` path at the top of `analyze.py` to point at your local copy of the sorted data folder.
2. Run:
   ```
   python analyze.py
   ```
   This iterates all samples, detects Dnd1 spots, segments each cell (reusing `masks/cell.tif` if it already exists, unless `recompute_masks = True` is set), computes per-cell features, and writes the outputs listed below.
3. (Optional, prior to/alongside analysis) `curate_cell_masks.ijm` is a Fiji/ImageJ macro for manual, blinded correction of `masks/cell.tif`. Run it in Fiji; it prompts interactively for the root data folder and only ever displays the blinded `uid` prefix, never the full sample name.

## Outputs

| Output file | Produced by | Description | Figure panel |
|---|---|---|---|
| `output/contact_sheet.pdf` | `analyze.py` (main loop + `fig_contact_sheet.savefig`) | Per-sample composite figure, one row per cell: merged Vasa/membrane/Dnd1, Vasa with cell contour, membrane with membrane-ring contour, raw Dnd1, DoG-filtered Dnd1 with spot contours, and the distance-from-membrane map | TODO |
| `output/dnd_membrane_sum_distribution_spots_stats.txt` | `analyze.py` (statistics block) | Text report of the statistical comparison (t-test or Mann-Whitney U, chosen by normality) of `dnd_membrane_sum_distribution_spots` between conditions `tmd-n_d` and `n_d`: test type, statistic, p-value, group sizes, significance stars | TODO |
| `output/representative_uids.txt` | `representative_images.export_representative_uids` (called from `analyze.py`) | For each of the 10 `PLOTTED_FEATURES` and each condition, lists the blinded `uid`s of samples falling in the 40th-60th percentile of that feature — used to pick representative images | TODO |
| `output/spot_distances_kde.pdf` | `analyze.py` | KDE plot of Dnd1 spot distance from the membrane (µm, x-axis) by condition (colored by condition) | TODO |
| `output/cell_features.pdf` | `analyze.py` | Strip + point plot (median, SD error bars) of 10 per-cell features by condition: `vasa_mean_intensity`, `membrane_mean_intensity`, `cell_area`, `dnd_mean_intensity`, `spot_number`, `spots_area`, `spots_threshold`, `spot_distances_mean`, `dnd_membrane_sum_distribution_all`, `dnd_membrane_sum_distribution_spots` | TODO |
| `<sample_dir>/masks/cell.tif` | `analyze.py` (main loop) | Cell segmentation mask (uint8), written into the sample's own folder outside this repo; reused on subsequent runs unless `recompute_masks = True` | TODO |
| `<sample_dir>/masks/spots.tif` | `analyze.py` (main loop) | Dnd1 spot segmentation mask (uint8), written into the sample's own folder outside this repo | TODO |

## Key files

- `analyze.py` — main pipeline: parses sample folders, detects Dnd1 spots, segments cells, computes distance-from-membrane and membrane/cytoplasm intensity features, runs statistics, and generates all plots.
- `segmenter.py` — wraps `micro_sam` for automatic cell segmentation (`SegmentationModel`), including checkpoint loading and GPU/MPS memory checks.
- `filtering.py` — Dnd1 spot detection primitives: difference-of-Gaussians filtering, MAD-based thresholding, connected-component labeling.
- `plotlib.py` — plotting helpers: pseudocolor single-channel display, additive multi-channel composites, scale bars.
- `utils.py` — statistics (normality check, t-test/Mann-Whitney selection, significance stars), `.lsm` pixel-size reading, image normalization, and mask utilities (membrane-ring extraction, center-object selection).
- `representative_images.py` — selects per-condition, per-feature "typical" sample UIDs for figure selection.
- `curate_cell_masks.ijm` — Fiji macro for blinded manual review/correction of cell masks.
