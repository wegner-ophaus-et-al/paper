# PatA RNA Localization

Quantifies whether Pateamine A (PatA), a translation inhibitor, changes the localization of `nanos` and `tdrd7` mRNA relative to Vasa+ germ granules in zebrafish primordial germ cells at 10 hpf, compared to DMSO-treated controls. For each imaged cell, mRNA intensity (RNAscope) is measured inside granules, in a granule periphery shell, and in the surrounding cytoplasm, and compared between conditions.

## Input data

Raw confocal `.lsm` files are not included in this repo (they're in `data/`, which is gitignored) and are available from the authors on request.

The analysis (`analysis.py`) expects a directory of per-sample folders, e.g. `data/RNAscope_tdrd7_nanos/collected_data/`, each named `{UID}__{condition}__{name}` (`condition` is `ctrl` or `exp`, translated to "DMSO Control" / "PatA Treated"), containing:

- `raw/granules.tif`, `raw/nanos-rna.tif`, `raw/tdrd7-rna.tif` — single-channel confocal images
- `segmentation/cell.tif`, `segmentation/granules.tif`, `segmentation/nanos_rna.tif`, `segmentation/tdrd7a_rna.tif` — instance segmentation masks

**Segmentation is not reproducible from this repo as-is.** `segmentations/segmenter.py` wraps `micro_sam` automatic instance segmentation, but it loads a trained model checkpoint from a path outside this repo. `proto/proto_segmentation.ipynb` documents the prototype approach used to generate the `raw/` and `segmentation/` folders, but most of its write calls are commented out and it isn't a clean runnable script.

`preprocessing/rearrange_imaging_repeat.py` organizes raw `.lsm` files into the `{UID}__{condition}__{name}` folder layout above (condition inferred from filename: `F125`/`dmso` → `ctrl`, else `exp`).

## Environment

No `requirements.txt`, `environment.yml`, or `pyproject.toml` is present in this repo. Based on the imports actually used:

- `marimo` (analysis.py is a marimo notebook/app)
- `numpy`, `pandas`, `scipy`
- `matplotlib`, `seaborn`
- `scikit-image` (`skimage`)
- `tifffile`
- `micro_sam` (only needed for `segmentations/segmenter.py`, i.e. the non-reproducible segmentation step)

You'll need to install these yourself; versions aren't pinned anywhere in the repo.

## How to run

1. Install the packages listed above.
2. Arrange your data under the directory layout described in "Input data" above.
3. In `analysis.py`, update the hardcoded `experiment_root` path (currently the author's local machine path) to point at your `collected_data` directory.
4. From the repo root, run the analysis:
   - `marimo edit analysis.py` (interactive), or
   - `python analysis.py` (headless)

## Outputs

| Output file | Produced by | Description | Figure panel |
|---|---|---|---|
| `{experiment_root}/summary_analysis.pdf` | `analysis.py`, final cell of the per-sample loop (`fig.savefig(...)`) | QC grid, 6 panels × N samples: raw granule/nanos-RNA/tdrd7-RNA channels and their segmentation overlays, one row per sample. Written into the data directory, not the repo. | TODO |
| `{sample_dir}/output/summary_images.png` | `analysis.py`, `process_sample()`, only when called with `save_images=True` (the default driver call in `analysis.py` uses `save_images=False`, so this is not produced by a default run) | Per-sample 4-panel figure: granules channel with periphery mask overlay, nanos RNA channel, tdrd7 RNA channel, and a merged RNA image, each with segmentation contours. | TODO |
| `quantifiaction.pdf` | `analysis.py`, `fig2.savefig("quantifiaction.pdf")` | 2×4 grid of violin/swarm plots comparing DMSO vs PatA, split by RNA type, for granule/periphery/cytoplasm intensity ratios, RNA–granule overlap fraction, RNA puncta count, mean granule size, and granule count. | TODO |
| `proto/output/quantification_paper.pdf` | `analysis.py`, `fig_paper.savefig(...)` | Single-panel violin/swarm plot of the cytoplasm/granule mean-intensity ratio by RNA type and condition, sized for the paper figure. | TODO |
| `proto/output/paper_figure_newstyle.pdf` | `analysis.py`, `fig_paper_new_style.savefig(...)` | Strip + point plot ("mRNA relocalization - PatA") of the cytoplasm/granule mean-intensity ratio by RNA type and condition, median with SD error bars, transparent background. | TODO |
| `proto/output/statistical_testing_results.csv` | `analysis.py`, `df_significance.to_csv(...)` | One row per (feature, RNA type) pair. Columns: `feature`, `rna_type`, `test_type` (t-test or Mann-Whitney U, chosen by normality test), `t_statistic`, `p_value`, `significance_symbol`. Compares DMSO Control vs PatA Treated. | TODO |

Note: `proto/quantifiaction.pdf` is also committed but no script in the repo currently writes to that exact path — it appears to be a leftover copy from an earlier run of the `analysis.py`/`proto/analysis.ipynb` cell (`fig2.savefig("quantifiaction.pdf")`) executed with a different working directory.

## Key files

- `analysis.py` — main analysis, a marimo notebook/app. Computes per-sample granule/periphery/cytoplasm masks (`generate_granule_periphery_mask`), extracts intensity and morphology features (`process_sample`), aggregates into a dataframe, computes intensity ratios, plots comparisons, and runs significance testing.
- `utils.py` — standalone helpers (mask computation, confocal pixel size extraction from `.lsm` metadata, statistical testing) — not currently imported by `analysis.py`, which defines its own local copies of some of the same logic.
- `colormaps.py` — custom black-to-color matplotlib colormaps (cyan, magenta, yellow, green, red) used for channel overlays.
- `segmentations/segmenter.py` — `SegmentationModel`/`SegmentationInstance` wrapper classes around `micro_sam` automatic instance segmentation.
- `preprocessing/rearrange_imaging_repeat.py` — sorts raw `.lsm` files into per-sample folders.
- `proto/` — earlier prototype versions of the analysis and segmentation notebooks, kept for reference; not the primary entry point.
