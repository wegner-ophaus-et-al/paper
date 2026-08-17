# Stress Granule Markers — Segmentation & Analysis

Quantifies stress granule formation in zebrafish embryos (10 hpf) expressing one of two fluorescent reporter constructs, F617 (eIF2s1a) or F618 (G3bp1), comparing control vs. heat-shock conditions. Per-cell images are segmented into granule/cytoplasm compartments with a fine-tuned Segment Anything (SAM) model, and per-cell area and intensity metrics are extracted and statistically compared between conditions.

## Pipeline

Manual (ImageJ/Fiji) and scripted (Python) steps alternate, in this order:

1. **`imageJ_macros/20260422_crop_masks.ijm`** (or `_2.ijm`) — opens each raw multichannel z-stack, lets the user manually select and mask individual cells from the max-intensity projection, and saves each selection as a mask in the sample's `masks/` folder.
2. **`crop_original_raw.py`** — for each per-cell mask, crops the raw stack to that mask's bounding box, hashes the cropped image to a short `uid`, and writes the cropped stack plus per-channel max-projections (`nls.tif`, `stress.tif`, `gra.tif`) to a new per-cell sample folder.
3. **`imageJ_macros/generate_cell_segmentations.ijm`** — opens the cropped `stress.tif`/`nls.tif` images and lets the user manually draw a whole-cell outline, saved as `masks/cell.tif`.
4. **`perform_granule_semgentations.py`** (main entry point) — for every sample: loads the three channel images and the cell mask, runs the fine-tuned SAM model (via `segmenter.py`/`processor.py`) to segment granules, computes areas and intensities, and generates the summary contact sheet, results table, plot, and statistics (see Outputs below).
5. **`summarize_control_vs_heatshock.py`** — reads the results CSV from step 4 and runs the same control-vs-heatshock statistical comparison across *every* numeric feature (not just the subset covered in step 4), plus per-condition cell counts.
6. **`imageJ_macros/adjust_bc_to_percentile.ijm`** (optional) — background-subtracts and contrast-adjusts a hand-picked set of representative samples to produce publication-style per-channel and merged PNGs.

## Input Data

Raw multichannel fluorescence z-stacks (channels: NLS nuclear reference, stress marker, granulito/granule marker), organized as `{construct_id}/{sample_name}/original_raw/*.tif`, one folder per imaged field of view. Manual per-cell masks (step 1 above) and cell outlines (step 3 above) are stored alongside the images as described in Pipeline.

**Data availability: TBD.** The scripts in this repo currently point at local/lab-internal paths (see Environment below); these need to be updated once the raw data has a permanent public location.

## Environment

**No environment or requirements file is included in this repo** — dependencies must be inferred from the imports below and installed manually:

- `numpy`, `scipy`, `pandas`
- `matplotlib`, `seaborn`
- `tifffile`, `scikit-image`, `Pillow`
- `tqdm`
- [`micro_sam`](https://github.com/computational-cell-analytics/micro-sam) (and its `torch` dependency) — used for SAM-based granule segmentation
- `psutil` (optional; used for MPS/Apple Silicon memory checks in `segmenter.py`, skipped if unavailable)

## How to Run

All data-root paths (e.g. in `crop_original_raw.py`, `perform_granule_semgentations.py`, `summarize_control_vs_heatshock.py`, and the top of each `.ijm` macro) are currently hardcoded local paths and must be edited to point at your own data before running. There is no CLI/config layer.

```bash
# 1 & 3: run the ImageJ macros interactively in Fiji (manual masking steps)

# 2: crop raw stacks to per-cell images
python crop_original_raw.py

# 4: main segmentation + analysis pipeline
python perform_granule_semgentations.py

# 5: extended cross-feature statistics
python summarize_control_vs_heatshock.py

# 6: optional representative-image export, run interactively in Fiji
```

`perform_granule_semgentations.py` requires a path to a fine-tuned `micro_sam` checkpoint (`SegmentationModel(path, model_type="vit_l_lm")`); this checkpoint is not included in the repo and is currently available on request.

## Outputs

`{data_root}` = the root data folder used by `perform_granule_semgentations.py` (contains one subfolder per construct, `F617`/`F618`). `{raw_data_root}` = the root used by `crop_original_raw.py`. `{figures_root}` = the folder used by `summarize_control_vs_heatshock.py`.

| Output file | Produced by | Description | Figure panel |
|---|---|---|---|
| `{data_root}/granule_segmentation_summary.pdf` | `perform_granule_semgentations.py` (main loop) | Contact sheet, one row per sample: granulito channel and stress-marker channel (both with granule-mask contour overlay), and NLS channel (with cell-mask contour overlay); each row labeled with sample `uid` and condition | TODO |
| `{data_root}/granule_segmentation_results.csv` | `perform_granule_semgentations.py` | One row per segmented cell. Columns: `uid`, `sample_name`, `construct_id`, `condition`, `total_cell_area`, `total_granules_area`, `total_cytoplasm_area`, `mean_int_granulito_{cell,granules,cytoplasm}`, `mean_int_stress_{cell,granules,cytoplasm}`, `mean_int_ratio_stress_cytoplasm_to_granule` (areas in pixels) | TODO |
| `{sample_folder}/masks/granules.tif` | `processor.py: process_sample()` | Per-sample granule instance segmentation mask (uint8) from the SAM model; reused instead of regenerated if already present, unless `resegmentation=True` | TODO |
| `{sample_folder}/visualization.png` | `processor.py: process_sample()` | Per-sample 3-panel QC figure (granulito, stress marker, NLS with mask contours); only written when `process_sample()` is called with `ax=None`, which does not happen in the current `perform_granule_semgentations.py` main loop (it always passes an axis) | TODO |
| `{data_root}/cytoplasm_granule_stress_pointplot.pdf` | `perform_granule_semgentations.py` | Strip + point plot of `mean_int_ratio_stress_cytoplasm_to_granule` (y) vs. `construct_id` (x: F618, F617), split by `condition` (control/heatshock); points show the median with SD error bars | TODO |
| `{data_root}/representative_ids_{construct_id}_{feature}.txt` | `perform_granule_semgentations.py` | Text list of sample `uid`s, per condition, whose value of `feature` falls between the 40th and 60th percentile for that construct/condition (i.e. representative, near-median samples) | TODO |
| `{data_root}/granule_segmentation_statistical_results.csv` | `perform_granule_semgentations.py` | One row per `construct_id` × `feature`. Columns: `Construct`, `Feature`, `P-value`, `Test type` (t-test or Mann-Whitney U, chosen by normality test), `significance` (stars), `Mean Control`, `Mean Heatshock`, `Percentage Change`, `Count Control`, `Count Heatshock` | TODO |
| `{raw_data_root}/cropped_images/{construct}/{uid}__{sample_name}-{idx}/images/raw_crop.tif` | `crop_original_raw.py` | Cropped multichannel raw z-stack for one cell, bounding-box-cropped from the original field of view using the manual selection mask | TODO |
| `{raw_data_root}/cropped_images/{construct}/{uid}__{sample_name}-{idx}/images/{nls,stress,gra}.tif` | `crop_original_raw.py` | Per-channel max-intensity projection of the cropped stack, normalized to 16-bit; `nls` = nuclear reference, `stress` = stress marker, `gra` = granulito (granule marker) | TODO |
| `{raw_data_root}/cropped_images/{construct}/{uid}__{sample_name}-{idx}/masks/{mask_filename}` | `crop_original_raw.py` | Copy of the manual per-cell selection mask (from step 1) used to define the crop | TODO |
| `{raw_data_root}/cropped_images/{construct}/{uid}__{sample_name}-{idx}/original_raw/{filename}` | `crop_original_raw.py` | Copy of the original uncropped raw stack for that field of view | TODO |
| `{sample_folder}/masks/cell.tif` | `imageJ_macros/generate_cell_segmentations.ijm` | Manually drawn whole-cell segmentation mask, drawn on the stress-marker/NLS composite | TODO |
| `{construct_folder}/{sample_name}masks/{original_filename}%{q}.tif` | `imageJ_macros/20260422_crop_masks.ijm` (or `_2.ijm`) | Manually drawn per-cell selection mask(s) used to crop individual cells out of a multi-cell field of view; `q` increments for each additional cell accepted from the same field | TODO |
| `{figures_root}/cell_counts_per_construct_condition.csv` | `summarize_control_vs_heatshock.py` | Cell counts per `construct_id` × `condition`, plus a `total` column; one row per construct | TODO |
| `{figures_root}/granule_segmentation_statistical_results_all_features.csv` | `summarize_control_vs_heatshock.py` | Same statistical comparison as `granule_segmentation_statistical_results.csv`, but run over every numeric feature in the results table (not just the three covered in the main pipeline) | TODO |
| `{export_dir}/{drug}/{nls,granulito,stress_protein}.png` | `imageJ_macros/adjust_bc_to_percentile.ijm` | Background-subtracted, percentile-contrast-adjusted single-channel PNG with scale bar; only produced when this macro is run manually on hand-picked representative samples | TODO |
| `{export_dir}/{drug}/merge.png` | `imageJ_macros/adjust_bc_to_percentile.ijm` | Composite merge PNG of the three adjusted channels | TODO |

`{construct_id}` / `{construct}` ∈ `F617`, `F618`. `{feature}` ∈ `mean_int_ratio_stress_cytoplasm_to_granule`, `total_cell_area`, `total_granules_area`.

## Key Files

| File | Role |
|---|---|
| `segmenter.py` | `SegmentationModel` (loads/runs the fine-tuned SAM model) and `SegmentationInstance` (single-image wrapper, used for model comparison during development, not in the main pipeline) |
| `processor.py` | `process_sample()` — per-sample segmentation, metrics, and QC panel |
| `perform_granule_semgentations.py` | Main entry point — orchestrates segmentation, results table, plots, and statistics |
| `crop_original_raw.py` | Crops raw multichannel stacks to per-cell images using the manual selection masks |
| `utils/utils.py` | Image I/O, normalization, upsampling, statistical tests, representative-sample selection |
| `utils_unused.py` | Earlier version of some `utils/utils.py` functions; not imported anywhere in the pipeline |
| `imageJ_macros/` | ImageJ/Fiji macros for the manual masking and representative-image steps |
