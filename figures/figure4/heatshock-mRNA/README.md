# Heat shock RNA localization analysis

Quantifies how heat shock affects the subcellular localization of maternal mRNAs (β-actin, *nanos3*) relative to germ granules (marked by Vasa) in zebrafish primordial germ cells (PGCs) at 12 hpf (Tg24 line), compared to somatic cells and to unstressed controls. Cells are segmented from confocal images, mRNA/granule channel intensities are measured per subcellular compartment (nucleus, cytoplasm, granules, whole cell), and cytoplasm-to-granule and cytoplasm-to-total intensity ratios are compared between conditions with violin/point plots and Mann-Whitney U / t-tests.

## Input data

Raw data are confocal `.lsm` stacks (RNAscope for β-actin/*nanos3* mRNA plus a Vasa/granule and a membrane or nucleus channel), one per imaged cell/embryo sample, organized on disk as:

```
<data_root>/<repeat_dir>/<uid>__<condition>-<sample_number>[__<date>]/
    original_raw/<condition>-<sample_number>.lsm
    raw/{nucleus,granule,bact_rna,nos_rna}.tiff        # channels split out of the .lsm
    segmentations/{nucleus,granule,cell}.tiff          # instance segmentation masks
    granule_coords.txt                                 # coordinates of the representative granule (from manual QC)
```

`condition` is one of `control`, `heatshock`, `somaControl`, `somaHeatshock` (the `soma*` samples are somatic-cell controls; the analysis script skips any sample whose condition contains `soma`).

**Data availability:** the raw imaging data is not included in this repository and has not yet been publicly deposited. The notebooks currently point at local/lab storage volumes (see "Key files" below) — update those paths to wherever the data is deposited before running.

## How to run it

No `requirements.txt` or `environment.yml` is included in this repository. From the imports used, you'll need at least: `marimo`, `numpy`, `pandas`, `matplotlib`, `seaborn`, `scipy`, `tifffile`, `tqdm`, `scikit-image`, and [`micro-sam`](https://github.com/computational-cell-analytics/micro-sam) (for segmentation only).

1. **Split raw stacks and run automatic segmentation** — `segment_20260114_images.ipynb` (repeat 1) / `segment_20260312_heatshock.ipynb` (repeats 2–3). Splits each multi-channel `.lsm` into single-channel TIFFs under `raw/`, then runs `SegmentationModel` (`segmentations/segmenter.py`, backed by `micro_sam.automatic_segmentation`) to produce initial instance segmentations under `segmentations/`. Edit the `root` path and the model checkpoint paths at the top of the notebook before running.
2. **Manually correct segmentations** — `20260116_screen-segmentation_heatshock.ijm` (repeat 1) / `20260312_screen-segmentation_heatshock.ijm` (repeats 2–3), run in Fiji/ImageJ. Steps through each sample, lets you remove spurious granule detections, mark a representative granule (saved to `granule_coords.txt`), and hand-draw nucleus/cell outlines, overwriting the `segmentations/*.tiff` files. Edit the hardcoded `root` path at the top of the macro, and the starting loop index if resuming a partially completed run.
3. **Measure intensities and compare conditions** — `analysis.py`, a [marimo](https://marimo.io) notebook (the current/canonical version of this analysis; `process_20260114_heatshock.ipynb` and `process_20260317_all_repeats.ipynb` are earlier, now-superseded versions of the same processing). Run with:
   ```
   marimo edit analysis.py
   ```
   Edit the `root` path in the second-to-last cell block to point at the sample data (a directory containing one subdirectory per repeat). It uses `utils.py` (`get_confocal_pixel_size`, `compute_masks`, `statistical_analysis`, `get_stat_stars`) to read pixel size from the original `.lsm` metadata, build cytoplasm/nucleus/granule/cell masks from the segmentations (cytoplasm = cell minus nucleus and granules, with a safety margin), measure mean/total marker intensity per compartment, and run condition comparisons.

## Outputs

| Output file | Produced by | Description | Figure panel |
|---|---|---|---|
| `raw/{channel}.tiff` | `segment_20260114_images.ipynb` / `segment_20260312_heatshock.ipynb` | Single-channel TIFF split out of the sample's raw `.lsm` stack, written into the external data directory (not this repo). `{channel}` ∈ `nucleus, granule, bact_rna, nos_rna` (repeat 1) or `granules` only (repeats 2–3, other channels currently commented out) | TODO |
| `segmentations/nucleus.tiff` | `segment_20260114_images.ipynb`, `SegmentationModel.segment` | Initial automatic instance segmentation mask of the nucleus channel, written into the external data directory | TODO |
| `segmentations/granules.tiff` | `segment_20260114_images.ipynb`, `SegmentationModel.segment` | Initial automatic instance segmentation mask of the granule channel, written into the external data directory | TODO |
| `granule_coords.txt` | `20260116_..._heatshock.ijm` / `20260312_..._heatshock.ijm` | Manual QC output: X/Y pixel coordinates (tab-separated, header `X\tY`) of the user-selected representative granule for the sample, written into the external data directory | TODO |
| `segmentations/granule.tiff` | `20260116_..._heatshock.ijm` / `20260312_..._heatshock.ijm` | Manually cleaned granule segmentation mask (spurious detections removed), overwritten into the external data directory | TODO |
| `segmentations/nucleus.tiff` | `20260116_..._heatshock.ijm` / `20260312_..._heatshock.ijm` | Manually drawn/refined nucleus mask, overwritten into the external data directory | TODO |
| `segmentations/cell.tiff` | `20260116_..._heatshock.ijm` / `20260312_..._heatshock.ijm` | Manually drawn whole-cell outline mask, written into the external data directory | TODO |
| `cytoplasm_to_granule_ratio_violinplot.pdf` | `analysis.py` | Split violin + swarm plot of cytoplasm-to-granule mean intensity ratio (unitless) per marker (`bactin_rna`, `nanos_rna`), split by condition (`control` vs `heatshock`) | TODO |
| `cytoplasm_to_granule_ratio_pointplot.pdf` | `analysis.py` | Strip + median point plot (alternate styling for figure assembly) of the same cytoplasm-to-granule ratio data as above | TODO |
| `cytoplasm_to_total_ratio_violinplot.pdf` | `analysis.py` | Split violin + swarm plot of cytoplasm intensity as a fraction of cytoplasm+granule intensity, per marker and condition | TODO |

Statistical test results (test type, statistic, p-value, significance stars, sample sizes) are printed to stdout by `analysis.py`, not written to a file.

`process_20260114_heatshock.ipynb` and `process_20260317_all_repeats.ipynb` are superseded versions of the `analysis.py` processing and write their own (now-unused) set of PDFs into the external data directory (e.g. `pgc_cytoplasm_granule_ratio_heatshock_effect.pdf`, `soma_cytoplasm_heatshock_effect.pdf`, and similar PGC/soma comparison plots), plus an optional per-sample `{sample_name}_overview.png` QC figure when `save_figures = True` is set in `process_20260114_heatshock.ipynb`.

## Key files

- `analysis.py` — canonical marimo analysis notebook (see above).
- `utils.py` — shared helpers: mask computation, confocal pixel size extraction, statistical tests, significance-star formatting.
- `segmentations/segmenter.py` — `SegmentationModel` / `SegmentationInstance` wrappers around `micro-sam` for automatic instance segmentation, including an upsample-segment-downsample path for small structures.
- `segment_20260114_images.ipynb`, `segment_20260312_heatshock.ipynb` — per-repeat raw-stack splitting and automatic segmentation.
- `20260116_screen-segmentation_heatshock.ijm`, `20260312_screen-segmentation_heatshock.ijm` — per-repeat interactive Fiji macros for manual segmentation QC/correction.
- `process_20260114_heatshock.ipynb`, `process_20260317_all_repeats.ipynb` — superseded exploratory versions of the `analysis.py` pipeline, kept for reference.
