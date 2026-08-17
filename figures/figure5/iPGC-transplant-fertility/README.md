# iPGC transplantation fertility replotting

Zebrafish were transplanted at the blastula stage with induced primordial germ cells (iPGCs) — either a full donor mix or donor iPGCs lacking *tdrd7a* — and outcrossed as adults to score their fertility and germline transmission of a dominant red-eye marker. This analysis takes the per-fish scoring data, computes fertilization rate, germline (marker) transmission rate, fertile/marker-positive proportions, and sex ratio per condition, and produces a single multi-panel summary figure. See `methods.md` for the full experimental methods and figure legend.

## Input data

`data/iPGC_transplantation_fertility_data.csv` — raw, manually compiled per-fish scoring data, one row per fish:

- `condition` — `wt`, `full_mix`, or `no_tdrd7a`
- `fish_id` — fish identifier (unique within a condition/stock, not globally)
- `stock_number` — donor/recipient stock (e.g. `H211`, `H213`)
- `fertilized`, `not_fertilized` — counts of fertilized vs. unfertilized eggs from the fish's clutch
- `red_eye_positive`, `red_eye_negative` — counts of F1 progeny with/without the dominant red-eye transgene (blank for `wt`, which carries no marker)

## How to run it

No environment file (`requirements.txt`/`environment.yml`/`pyproject.toml`) is provided. The script requires Python 3 with `pandas`, `matplotlib`, and `seaborn` installed, and a system installation of the `Arial` font (set via `mpl.rcParams["font.family"]`).

```
pip install pandas matplotlib seaborn
python replotting.py
```

## Outputs

| Output file | Produced by | Description | Figure panel |
|---|---|---|---|
| `data/iPGC_transplantation_fertility_data_processed.csv` | `replotting.py` | Per-fish data with computed columns added: `percent_fertilized`, `percent_dominant_marker`, `is_fertile`, `is_dominant_marker`. Fish with partial (not 0% or 100%) dominant-marker transmission are excluded. | TODO |
| `data/iPGC_transplantation_fertility_data_summary.csv` | `replotting.py` (`df.groupby("condition").describe()`) | Summary statistics (count, mean, std, quartiles, etc.) of the processed columns, grouped by `condition`. | TODO |
| `figures/iPGC_transplantation_fertility_replot.pdf` | `replotting.py` | 5-panel summary figure: (1) strip/point plot of fertilization rate (%) by condition; (2) strip/point plot of germline (dominant marker) transmission rate (%) by condition; (3) bar plot of proportion fertile per condition with 95% CI; (4) bar plot of proportion dominant-marker-positive per condition with 95% CI; (5) bar plot of sex ratio (% males) per condition. | TODO |

## Key files

- `replotting.py` — loads the raw CSV, computes derived metrics, writes the processed/summary CSVs, and generates the figure.
- `methods.md` — methods text and figure legend describing the experiment and this figure as they appear in the manuscript.
