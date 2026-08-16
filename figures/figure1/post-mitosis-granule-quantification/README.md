# Post-mitosis germ granule quantification

Quantifies germ granules (RNP granules, visualized with a fluorescent reporter) in zebrafish germ cells (~10 hpf) in the minutes following cell division. For each imaged cell, granule area is tracked frame-by-frame (5 min intervals) starting at `time = 0`, defined as the first frame after completion of cytokinesis. The analysis compares wild-type cells against `tdrd7a` knockdown (`tdrd7aMO`) cells, computing granule count, area, and volume (estimated from area assuming spherical granules) over time, and testing whether granule count at a given timepoint differs between conditions.

## Input data

- `data/all_data.xlsx` — the actual input read by the notebook. One row per cell per timepoint, with columns `condition` (`wt`/`tdrd7a`), `stage`, `date`, `event-number`, `cell-id`, `time` (minutes relative to cytokinesis), and a variable number of `area1`, `area2`, … columns holding the measured area (µm²) of each individual granule segmented in that cell/timepoint.
- `data/wt.xlsx`, `data/tdrd7aMO.xlsx` — earlier per-condition versions of the same data, superseded by `data/all_data.xlsx`. Not read by any script in this repo.
- `Quantification movies mitosis_new.xlsx` — raw per-movie granule area measurements (one sheet per imaged movie, plus summary sheets). This is the upstream source for `data/all_data.xlsx`, but the compilation step is not scripted in this repo.

## How to run

No environment/requirements file is included in this repo. The notebook depends on `pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`, and `openpyxl` (for reading `.xlsx`), run in Jupyter.

```
jupyter notebook 20240610_data-prep_plotting.ipynb
```

Run all cells top to bottom. The notebook reads `data/all_data.xlsx`, derives per-cell/per-timepoint granule metrics, and writes the figures below into `figures/`.

## Outputs

| Output file | Produced by | Description | Figure panel |
|---|---|---|---|
| `figures/granules_post-mitosis.pdf` | `20240610_data-prep_plotting.ipynb`, cell with `plt.savefig('figures/granules_post-mitosis.pdf', ...)` | 3-panel figure (wt vs. tdrd7a, mean ± SD): granule count vs. time after mitosis (min); relative total granule volume (volume at t / volume at t=0) vs. time; mean granule area (µm²) vs. time. X-axis restricted to 0–60 min. | TODO |
| `figures/20260607_granules_count_post_mitosis_for_fig1.pdf` | `20240610_data-prep_plotting.ipynb`, cell with `plt.savefig('figures/20260607_granules_count_post_mitosis_for_fig1.pdf', ...)` | Single-panel figure: granule count vs. time after mitosis (min), 0–60 min, wt vs. tdrd7a, mean ± SD. | TODO |

A Mann-Whitney U test comparing granule count between `wt` and `tdrd7a` at `time = 60` is printed inline in the notebook but not written to disk.

## Key files

- `20240610_data-prep_plotting.ipynb` — loads `data/all_data.xlsx`, converts per-granule areas to volumes (assuming spherical granules), computes per-cell/timepoint summary metrics (mean/median/max/min/sum area and volume, granule count, granule density, relative total volume), generates the figures above, and runs the granule-count Mann-Whitney U test.
