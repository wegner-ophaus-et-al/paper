# CHX germ granule mRNA localization analysis

Analysis of zebrafish primordial germ cells (10 hpf), comparing cycloheximide
(CHX) treatment against DMSO control. Cells are immunostained for Vasa
(germ granule marker) with RNAscope probes for *nanos* and *tdrd7a* mRNA.
The pipeline segments germ granules, cell bodies, and RNA puncta, then
quantifies how *nanos*/*tdrd7a* mRNA partitions between granules, a
granule-periphery shell, and the surrounding cytoplasm, to test whether
blocking translation (CHX) shifts mRNA localization relative to DMSO
controls.

## Input data

Each sample is a single-cell confocal image (`.lsm`, 3 channels: Vasa,
nanos RNAscope, tdrd7a RNAscope) acquired at 10 hpf under either CHX or
DMSO treatment. Raw imaging data is not included in this repository (not
yet publicly deposited) and the scripts currently point at a local path
(`/Users/julian/local_files/chx_kim`) that will need to be changed to
wherever the data actually lives.

Expected on-disk layout per sample, built up by the scripts below:

```
<sample_uid>__<condition>__<name>/
├── original_raw/<name>.lsm        # raw 3-channel confocal stack
├── imgs/{vasa,nanos,tdrd7a}.tif    # extracted single-channel images
├── segmentations/
│   ├── vasa.tif                   # granule segmentation (labeled)
│   ├── nanos.tif                  # nanos RNAscope puncta segmentation
│   ├── tdrd7a.tif                 # tdrd7a RNAscope puncta segmentation
│   └── cell.tif                   # whole-cell mask
└── output/summary_images.pdf
```

`condition` is `chx` or `ctrl` (mapped to "CHX Treated" / "DMSO Control").
`segmentations/cell.tif` (whole-cell outline) is **not** produced by any
script in this repo — it must be supplied or created by a separate,
manual step before running the ImageJ refinement macro or the analysis
script.

## How to run it

**Environment**: no `requirements.txt`/`environment.yml` is provided.
The scripts import `numpy`, `scipy`, `scikit-image`, `matplotlib`,
`seaborn`, `pandas`, `tifffile`, and (for segmentation only) `gamgee`
(`gamgee.segmenter.SegmentationModel`/`SegmentationInstance`), a
segmentation package/model that will be published as a companion repo
and is not included here. Install these yourself; there is no pinned
version list.

Steps, in order:

1. **`20260302_data_management.py`** — walks a directory of raw
   acquisition days/repeats, finds `.lsm` files, and copies them into a
   flat `data/` folder, renaming each to
   `{condition}__{day}_{repeat}_{filename}.lsm` (`condition` is `chx` if
   `"chx"` appears in the parent folder name, otherwise `ctrl`).
2. **`20260302_segmentation.py`** — for each sample under `data/`, loads
   the `.lsm` stack, runs the `gamgee` SAM-based segmentation model on
   each of the three channels (Vasa, nanos, tdrd7a), and writes the raw
   channel images to `imgs/` and the resulting segmentations to
   `segmentations/`. Also saves a combined contact sheet.
3. Supply/create `segmentations/cell.tif` (whole-cell mask) for each
   sample — not automated in this repo.
4. **`20260302_screen_segmentation.ijm`** (Fiji/ImageJ macro) — opens
   each sample's Vasa image and granule segmentation side by side for
   manual review, lets you erase unwanted granule regions, and lets you
   manually refine `segmentations/cell.tif`, overwriting both files.
5. **`20260303_analysis.py`** — for each sample, builds granule/periphery/
   cytoplasm masks from the segmentations, computes per-cell intensity
   and morphology features, aggregates across samples, runs statistical
   comparisons between CHX and DMSO, and writes the tables/plots listed
   below.

## Outputs

All paths are relative to the repo root and land in `output/`, which is
gitignored.

| Output file | Produced by | Description | Figure panel |
|---|---|---|---|
| `output/quantifiaction.pdf` | `20260303_analysis.py` (`fig2`) | 2×4 grid of violin/swarm plots per RNA type and condition: granule/cytoplasm, granule-center/cytoplasm, periphery/cytoplasm, and periphery/granule-center mean-intensity ratios; RNA–granule overlap fraction; RNA puncta count; mean granule size; granule count. | TODO |
| `output/quantification_paper.pdf` | `20260303_analysis.py` (`fig_paper`) | Single split-violin/swarm plot of the cytoplasm/granule mean-intensity ratio, by RNA type (nanos vs tdrd7a) and condition (CHX vs DMSO). | TODO |
| `output/quantification_paper_newstyle.pdf` | `20260303_analysis.py` (`fig_paper_newstyle`) | Strip + point plot (median, SD error bars) of the cytoplasm/granule mean-intensity ratio by RNA type and condition, in the paper's Arial/6pt styling. | TODO |
| `output/quantification_results.csv` | `20260303_analysis.py` (`df.to_csv`) | Per-sample-per-RNA-type row of all computed features (granule count/sizes, Vasa intensity, RNA intensities in granules/periphery/cytoplasm/granule-centers, puncta counts, overlap areas, and derived ratios). Key columns: `UID`, `Condition`, `RNAType`, `MeanRNA_Cytoplasm_Granules_Ratio`. | TODO |
| `output/statistical_testing_results.csv` | `20260303_analysis.py` (`get_signifcances` / `df_significance.to_csv`) | One row per (feature, RNA type): CHX vs DMSO comparison (`test_type` is t-test or Mann-Whitney U depending on normality), `t_statistic`, `p_value`, `significance_symbol`. | TODO |
| `output/representative_samples.txt` | `20260303_analysis.py` (`get_representative_ids`) | For each feature and condition, sample `UID`s falling within the 40th–60th percentile of that feature — candidates for representative images. | TODO |
| `{sample_dir}/output/summary_images.pdf` | `20260303_analysis.py` (`process_sample`, only if `save_images=True`) | Per-sample 4-panel figure: Vasa image with periphery mask and cell outline, nanos RNA image, tdrd7a RNA image, and a merged nanos+tdrd7a overlay with granule contour. Off by default in the current run (`save_images=False`). | TODO |
| `{sample_dir}/segmentations/{vasa,nanos,tdrd7a}.tif` | `20260302_segmentation.py` | Per-channel segmentation masks from the `gamgee` model, later refined for granules by the ImageJ macro. | TODO |
| `{sample_dir}/imgs/{vasa,nanos,tdrd7a}.tif` | `20260302_segmentation.py` | Per-channel raw images extracted from the `.lsm` stack. | TODO |
| `{sample_dir}/segmentations/cell.tif` | ImageJ macro `20260302_screen_segmentation.ijm` (refinement only) | Whole-cell mask, manually refined; initial mask is created outside this repo. | TODO |
| `data/{condition}__{day}_{repeat}_{filename}.lsm` | `20260302_data_management.py` | Raw stack copied and renamed from the acquisition folder tree. `condition` is `chx` or `ctrl`. | TODO |
| `{date}_contactsheet.pdf` (written to the working directory) | `20260302_segmentation.py` | Grid of Vasa/nanos/tdrd7a images with segmentation contours overlaid, one row per sample, for a quick QC pass over the whole segmented set. `{date}` is today's date (`YYYYMMDD`) at run time. | TODO |

## Key files

- `20260302_data_management.py` — collects raw `.lsm` files into `data/`.
- `20260302_segmentation.py` — runs `gamgee` segmentation, writes `imgs/` and `segmentations/`.
- `20260302_screen_segmentation.ijm` — Fiji macro for manual review/refinement of granule and cell segmentations.
- `20260303_analysis.py` — feature extraction, statistics, and figure generation; the main analysis entry point.
- `merge_utils.py` — `merge_images_with_cmaps`, used to overlay two single-channel images with independent colormaps.
- `colormaps.py` — black-to-color matplotlib colormaps (cyan/magenta/yellow/green/red) used for channel display.
