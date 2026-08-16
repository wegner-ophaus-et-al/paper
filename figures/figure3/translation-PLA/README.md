# PLA translation-site analysis (duolink-PLA)

Analysis of proximity ligation assay (PLA) images in zebrafish primordial germ cells, used to detect sites of active translation of *vasa* and β-actin RNA (nascent Vasa / β-actin protein close to the ribosome exit site). For each imaged cell, the pipeline segments the cell and its germ granules, detects individual PLA puncta ("spots"), and measures spot counts, spot intensity, and spot-to-granule distance, comparing the drug conditions DMSO (control), CHX (cycloheximide), and PatA (Pateamine A).

## Input data

Expected under `data/`:

- A raw-data folder of per-cell image sets (currently `data/2026-05-26 RNAscope after PLA/`, gitignored — **not committed to this repo; data availability is still TBD**). Each subfolder is named `{uid}__{rna_probe}-{condition}-{sample_id}` (e.g. `f245e072__vasa-PatA-6`) and contains:
  - `original_raw/{name}.lsm` — raw 3-channel confocal stack (channel 0 = granule, 1 = PLA, 2 = RNA probe)
  - `images/{granule,pla,vasa|actin}.tif` — per-channel exported TIFFs
  - `segmentations/{granule,cell}.tif` — cell and granule segmentation masks
- `data/representative_uids.txt` — UIDs of representative cells per condition, used to pick example images.

`Cell` objects (`cell.py`) are discovered automatically in `main.py` by globbing `data/**/*.lsm`, so any correctly named/structured folder under `data/` will be picked up.

## How to run

**No environment file is committed** (no `requirements.txt`, `environment.yml`, or `pyproject.toml`) — dependency versions are not pinned anywhere in the repo. Based on the imports, you'll need: `numpy`, `scipy`, `scikit-image`, `tifffile`, `pandas`, `seaborn`, `matplotlib`, `tqdm`.

The default pipeline also needs `filtering.py` (spot-detection core: `segment_spots`, `mad_threshold`, `dog_filter`), which is currently excluded via `.gitignore` — it needs to be added to the repo for `spots.py`/`cell.py` to import successfully.

`main.py` additionally imports a SAM-based granule-segmentation package, `gamgee` (`from gamgee.segmenter import SegmentationModel`, loaded via a hardcoded `sys.path` entry). This package will be published separately and is **not required for a default run** — segmentation only runs if `segment_new=True` is set in `main.py`; by default (`segment_new=False`) the pipeline reads precomputed segmentation masks from `segmentations/` instead.

Once dependencies are available, run from the repo root:

```
python main.py
```

This processes every cell folder found under `data/`, computes measurements, and writes all outputs listed below into `data/`.

Two ImageJ macros are used as separate, manual/interactive steps (not invoked by `main.py`):
- `20260610_screen_segmentation.ijm` — opens each cell's images for manual review/refinement of the cell segmentation mask.
- `adjust_bc_to_percentile.ijm` — background-subtracts and contrast-adjusts channels for representative figure images, exporting PNGs per condition; the `drug` variable at the top of the script must be changed manually per condition/variant to reprocess.

Both macros currently have absolute local file paths hardcoded at the top and need to be edited to your own paths before running.

## Outputs

All outputs are written by `main.py`.

| Output file | Produced by | Description | Figure panel |
|---|---|---|---|
| `data/contactsheet_segmentation_results.pdf` | `main.py` (via `Cell.plot_segmentations`) | Per-cell contact sheet of segmentation results (cell/granule masks overlaid on images) | TODO |
| `data/contactsheet_generated_masks.pdf` | `main.py` (via `Cell.plot_images` / `generate_various_masks`) | Per-cell contact sheet of derived masks (granule periphery, pseudo-nucleus, cytoplasm, nucleo-cytoplasm) | TODO |
| `data/contactsheet_pla_spots.pdf` | `main.py` (via `Cell.plot_spot_detection`) | Per-cell contact sheet of detected PLA spots overlaid on the PLA channel | TODO |
| `data/measurements_all_cell.csv` | `main.py` | Per-cell, per-segmentation, per-channel, per-statistic intensity measurements (sum/mean/median) for all imaged cells, unfiltered. One row = one (cell, segmentation region, channel, statistic) combination. | TODO |
| `data/measurements_bright_cells.csv` | `main.py` | Same as above, restricted to cells passing the brightness QC filter (mean granule intensity > `brightness_threshold`, default 9000). Only written if `brightness_threshold` is set. | TODO |
| `data/figures/{channel}_measurements.pdf` | `main.py` (`compare_conditions` / `utils.statistics`) | Per-channel intensity-ratio comparison plot across DMSO/CHX/PatA. `{channel}` ∈ `granule`, `pla`, `rna_type`. | TODO |
| `data/figures/pla_spots_summary_{rna_tpe}.pdf` | `main.py` (via `spots.plot_spot_summary`) | Spot count and mean spot intensity per cell, by condition, split by RNA probe. `{rna_tpe}` ∈ `actin`, `vasa`. | TODO |
| `data/spot_count_statistics.csv` | `main.py` (via `utils.compare_conditions`) | Summary statistics (DMSO vs. CHX, DMSO vs. PatA) for spot counts/intensity — one row per comparison/condition. | TODO |
| `data/figures/pla_spot_distance_to_granule_kde.pdf` | `main.py` (via `spots.plot_distance_kde`) | KDE of PLA spot distance to nearest granule (µm), DMSO condition, split by RNA probe. | TODO |
| `data/figures/pla_spot_distance_to_granule_vs_intensity_hex.pdf` | `main.py` (via `spots.plot_distance_vs_intensity_hex`) | Hexbin joint plot of spot distance to granule (µm) vs. spot intensity, DMSO condition. | TODO |
| `data/statistical_results.txt` | `main.py` (via `utils.statistics`) | Text summary of all condition-comparison statistical tests (normality test → t-test or Mann-Whitney U, with significance stars) run during the analysis. | TODO |

## Key files

- `main.py` — pipeline entry point: discovers cells, runs segmentation/measurement/spot-detection per cell, aggregates results, produces all plots and tables above.
- `cell.py` — `Cell` class: loads raw/processed images and segmentations for one cell, builds derived masks (granule periphery, pseudo-nucleus, cytoplasm), runs spot detection and per-region intensity measurements.
- `spots.py` — PLA spot detection and evaluation (`detect_and_evaluate`), and spot-summary/distance plotting functions.
- `filtering.py` — core spot-segmentation algorithm (difference-of-Gaussians filtering, MAD-based thresholding). **Currently not committed** (see "How to run" above).
- `utils.py` — shared helpers: pseudo-nucleus convex-hull mask, area-ratio calculation, statistical tests, and `.lsm` pixel-size reading.
- `20260610_screen_segmentation.ijm`, `adjust_bc_to_percentile.ijm` — manual ImageJ macros for segmentation refinement and figure image export (see "How to run").
