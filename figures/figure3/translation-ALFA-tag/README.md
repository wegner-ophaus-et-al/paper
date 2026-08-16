# ALFA-tag confocal analysis (PatA)

Quantifies ALFA-tag signal in confocal images of cells segmented into
nucleus, granule, and cytoplasm compartments, comparing a PatA-treated
condition against Control. For each cell it measures per-compartment
ALFA-tag intensity, detects individual ALFA-tag spots and their distance to
the nearest granule, and runs Control-vs-PatA statistics on all measured
features.

## Input data

The raw data is not yet public (deposition/accession TBD). The pipeline
expects a root directory (set as `root` in `analyze.py`, currently hardcoded
to a local path — edit before running) containing one folder per sample:

```
<root>/<uid>__<sample_name>/original_raw/<sample_name>.lsm
<root>/<uid>__<sample_name>/masks/cell.tif
<root>/<uid>__<sample_name>/masks/nucleus.tif
<root>/<uid>__<sample_name>/masks/granule.tif
```

- `uid` is an 8-hex-character sample ID; `sample_name` is free text and
  determines the condition — folder names containing `"pata"`
  (case-insensitive) are labelled `PatA`, everything else `Control`.
- Each `.lsm` is a 4-channel confocal stack: channel 0 = Vasa, 1 = membrane,
  2 = mCherry (target expression), 3 = ALFA-tag.
- The `masks/*.tif` instance segmentation masks are expected to already
  exist. `analyze.py` can compute them via `segmenter.py`
  (a `micro_sam`-based `SegmentationModel`), but in the current committed
  state that model loading is commented out (`sm = None`, `sm_g = None`),
  so as-is the script requires precomputed masks.

## Environment and how to run
1. (Optional) Curate/correct cell masks in Fiji/ImageJ with the
   Bio-Formats plugin installed:
   ```
2. Edit `root` in `analyze.py` to point at your local data directory, then
   run the main pipeline:
   ```
   python analyze.py
   ```
3. Select representative sample UIDs per feature/condition from the
   resulting `output/cell_features.csv`:
   ```
   python select_representative_uids.py
   ```
4. (Optional, manual) `adjust_bc_to_percentile.ijm` is a one-off Fiji macro
   for exporting brightness/contrast-adjusted PNGs of a single
   representative sample for figure panels. `export_dir`, `drug`, and the
   channel numbers are hardcoded and meant to be edited per sample before
   running.

## Outputs

| Output file | Produced by | Description | Figure panel |
|---|---|---|---|
| `output/contact_sheet.pdf` | `analyze.py` | Per-sample QC grid (one row per sample): Vasa with granule contour, membrane with cell contour, mCherry with nucleus contour, raw ALFA-tag, Vasa/mCherry merge, compartment mask contours, and detected ALFA-tag spots. | TODO |
| `output/cell_features.csv` | `analyze.py` | One row per cell. Key columns: `uid`, `condition`, `pixel_size`, per-compartment ALFA-tag intensity (`mean_alfa_cytoplasm`, `sum_alfa_granules`, `cv_alfa_granule_periphery`, etc.), spot count/area/intensity (`number_of_spots`, `area_of_spots`, `sum_intensity_alfa_spots`, `mean_intensity_alfa_spots`), `mean_distance_to_granules` (µm), and cytoplasm:granule intensity ratios (`MeanRNA_Cytoplasm_Granules_Ratio`, `SumRNA_Cytoplasm_Granules_Ratio`). | TODO |
| `output/stats_results.csv` | `analyze.py` | One row per measured feature: `test_type` (t-test or Mann-Whitney U, chosen by a normality test), `n_control`/`n_pata`, `statistic`, `p_value`, `stars`, for Control vs PatA. | TODO |
| `output/basic_features.pdf` | `analyze.py` | 3×3 grid of strip+point plots (Control vs PatA) for mean, sum, and CV of ALFA-tag intensity in the cytoplasm, granules, and granule periphery. | TODO |
| `output/spots_features.pdf` | `analyze.py` | Strip+point plots (Control vs PatA) for spot count, spot area, summed and mean spot intensity, mean distance to granules, nucleus target expression, and summed cytoplasm:granule ratio. | TODO |
| `output/kde_distance_to_granules.pdf` | `analyze.py` | KDE of per-cell mean distance to nearest granule (µm, x-axis), Control condition only. | TODO |
| `output/spot_intensity_vs_distance_to_granules.pdf` | `analyze.py` | Per-spot hexbin joint plot: distance to nearest granule (µm, x-axis) vs. spot ALFA-tag intensity (y-axis), Control condition only. | TODO |
| `output/representative_uids.txt` | `select_representative_uids.py` | For each feature and condition, lists sample UIDs whose value falls in the 40th–60th percentile of `cell_features.csv` — candidates for representative images. | TODO |
| `<root>/<uid>__<name>/masks/{cell,nucleus,granule}.tif` | `analyze.py`; `cell.tif` also editable via `curate_cell_masks.ijm` | Per-sample instance segmentation masks, written back into the input data tree (not under `output/`) when not already present or when `recompute_masks=True`. | TODO |
| `<export_dir>/<drug>/{granules,alfa,merge}.png` | `adjust_bc_to_percentile.ijm` | Brightness/contrast-adjusted per-channel and composite PNGs for one hand-picked representative sample. `export_dir` and `drug` are hardcoded in the macro and edited per use; output path is outside this repo. | TODO |

## Key files

- `analyze.py` — main pipeline: segmentation loading, compartment masks,
  spot detection, feature measurement, statistics, and all plots.
- `segmenter.py` — `SegmentationModel`/`SegmentationInstance` wrappers
  around `micro_sam` automatic instance segmentation.
- `filtering.py` — spot detection: difference-of-Gaussians filter
  (`dog_filter`) plus a robust MAD-based threshold (`segment_spots`);
  `sweep_k` is a threshold-calibration helper not currently called from
  `analyze.py`.
- `utils.py` — geometry helpers (ring masks, mask confinement, center-object
  selection), LSM pixel-size extraction, and the statistics helpers
  (`statistical_analysis`, `get_stars`).
- `plotlib.py` — single-color-channel and multi-channel merge image
  plotting with custom colormaps, plus a scale-bar helper.
- `select_representative_uids.py` — picks representative sample UIDs per
  feature/condition from `output/cell_features.csv`.
- `curate_cell_masks.ijm` / `adjust_bc_to_percentile.ijm` — Fiji/ImageJ
  macros for manual mask correction and figure-panel image export,
  respectively (see "How to run" above).
