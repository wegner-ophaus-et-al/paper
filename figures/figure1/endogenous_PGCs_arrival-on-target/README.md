# Tdrd7a MO ectopicity

Quantifies primordial germ cell (PGC) migration defects in zebrafish embryos after morpholino (MO) knockdown of `tdrd7a`, compared to a control morpholino. For each embryo, PGCs are scored as either correctly localized to the gonad ("on target") or found outside it ("ectopic"), and the analysis compares the resulting per-embryo ectopicity percentage and total PGC number between the `ctrl` and `kd` (tdrd7a MO) conditions, with statistical testing.

## Input data

- `Tdrd7a_ectopicity.csv` — the raw, manually scored data. One row per embryo, with columns `Date`, `Embryo Nr`, `Gonad` (count of PGCs at the gonad), `Ectopic` (count of PGCs outside the gonad), and `condition` (`ctrl` or `kd`).
- `Tdrd7a ectopicity.xlsx` — the original spreadsheet the CSV was derived from. It is excluded from version control (see `.gitignore`) and is not required to run the analysis.

## How to run it

No environment file (`requirements.txt` / `environment.yml` / `pyproject.toml`) is included in this repo, and package versions are not pinned. The script requires:

```
pip install pandas matplotlib seaborn scipy
```

Then, from the repo root:

```
python plot.py
```

This reads `Tdrd7a_ectopicity.csv`, computes per-embryo percentages, runs the statistical comparison, and writes the outputs below into `figure/`.

## Outputs

| Output file | Produced by | Description | Figure panel |
|---|---|---|---|
| `figure/ePGCs_ectopicity_number.pdf` | `plot.py` | Two-panel strip/point plot (median ± SD) comparing `ctrl` vs `kd` for (1) `on_target` — percent of PGCs found in the gonad per embryo (y-axis, 0–101%) and (2) `total_number` — total PGCs counted per embryo (Gonad + Ectopic) | TODO |
| `figure/stats.txt` | `plot.py` | Text summary of the `ctrl` vs `kd` statistical comparison for `on_target` and `total_number`: test used (t-test or Mann-Whitney U, chosen via `utils.statistical_analysis` based on normality), p-value, significance stars, and n per group | TODO |

## Key files

- `plot.py` — the only script that runs; loads the CSV, derives `total_number`, `ectopicity`, and `on_target` columns, generates the figure, and writes the stats summary.
- `utils.py` — shared helpers: `statistical_analysis` (picks t-test vs Mann-Whitney U by normality test) and `get_stars` (p-value to significance stars).
- `plotting.py` — a larger set of plotting functions (granule/PGC-marker profiles, fold-change bar plots, ridge plots, etc.) for markers/stages not present in this repo's data. It is not imported or used by `plot.py`; it is unused leftover code carried over from another analysis and does not apply to the ectopicity data here.
- `ePGCs_ectopicity_number.pdf` (repo root) — a stale duplicate from before `plot.py` started writing into `figure/`; not produced by the current code.
