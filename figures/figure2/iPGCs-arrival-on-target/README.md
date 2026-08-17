# iPGC Tdrd7a Ectopicity Analysis

Quantifies where primordial germ cells (PGCs) end up in injected embryos, comparing embryos raised **with** vs **without** Tdrd7a. For each embryo, cells are counted in three regions (`Gonad`, `Ectopic`, `HighLigandArea`); `Gonad` and `HighLigandArea` are combined into "on-target" arrival, which is contrasted with `Ectopic` (mislocalized) placement. The analysis computes per-embryo percentages, plots the distribution by condition, and runs a statistical comparison of ectopic PGC rates between conditions.

## Input data

`ectopicity_data.csv` — manually counted PGCs per embryo, one row per embryo:

| Column | Meaning |
|---|---|
| `Date` | Experiment date |
| `Condition` | `with Tdrd7a` or `without Tdrd7a` |
| `Gonad` | Cell count in the gonad (on-target) |
| `Ectopic` | Cell count outside the target region |
| `HighLigandArea` | Cell count in a secondary on-target region |

## How to run

Requires `pandas`, `seaborn`, `matplotlib`, `scipy`, and `marimo` (no `requirements.txt`/`environment.yml`/`pyproject.toml` is present in this repo, so versions are unpinned — install these with `pip` as needed).

```bash
marimo edit main.py   # interactive
# or
marimo run main.py    # headless, regenerates all outputs
```

Running `main.py` end-to-end reads `ectopicity_data.csv`, reshapes it, computes per-embryo percentages, generates both plots, runs the statistical test, and writes `stats.txt`.

## Outputs

| Output file | Produced by | Description | Figure panel |
|---|---|---|---|
| `ectopicity_violinplot.pdf` | `main.py`, violin plot cell | Split violin plot of `Percent` (y, % of counted cells in that embryo) by `Region` (x: `Gonad`/`Ectopic`/`HighLigandArea`), split by `Condition`. Only embryos with total count > 10 (`MIN_CELL_NUMBER`) are included. Not committed to git (gitignored). | TODO |
| `ectopicity_paper_plot.pdf` | `main.py`, final plotting cell | Strip + point plot of on-target `Percent` ("OnTarget" = `Gonad` + `HighLigandArea`) by `Condition`, titled "On target iPGCs at 24 hpf"; points show per-embryo values, markers show median ± SD. | TODO |
| `stats.txt` | `main.py`, stats cells | Text summary of the ectopic-PGC comparison between conditions: sample sizes, test used (t-test or Mann-Whitney U, chosen via `utils.statistical_analysis` based on normality), p-value, and significance stars. | TODO |

## Key files

- `main.py` — marimo notebook; the sole entry point. Loads and reshapes the data, computes percentages, generates both plots, and runs the statistical test.
- `utils.py` — `statistical_analysis()` (picks t-test or Mann-Whitney U based on a normality test) and `get_stars()` (p-value → significance stars).
