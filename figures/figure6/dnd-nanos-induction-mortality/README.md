# TMD-nanos induction mortality

Analysis of embryo survival/malformation counts from a zebrafish experiment inducing germ cell fate in primordial germ cells (PGCs) using a membrane-bound nanos3 construct (tmd-nanos). Embryos injected under four different RNA combinations were scored as dead, severely malformed, or normal, and this analysis computes and plots the fraction of severe/dead outcomes per condition and tests for differences between conditions.

## Input data

`tmd-nanos_conditions.csv` — the raw scored data (manually entered counts), one row per date/condition/category/repeat combination. Columns:

- `date` — date the repeat was scored.
- `condition` — injection condition:
  - `n_d` — induction with dead end (dnd1) + (non-membrane-bound) nanos3.
  - `tmd-n_d` — induction with membrane-bound nanos3 + free dnd1.
  - `tmd-n_stuffer` — membrane-bound nanos3 + a control (stuffer) RNA, no induction.
  - `wt` — wild type, injected with a membrane marker + control RNA.
- `category` — outcome: `dead`, `sev_malformed`, or `normal`. Embryos were assessed at 24 hpf and confirmed at 48 hpf; severely malformed embryos are pooled with dead ones in the analysis because all had died by 48 hpf.
- `count` — number of embryos in that date/condition/category/repeat bucket.
- `repeat` — experimental repeat number (1–3; each condition has 3 repeats).

## How to run

```bash
pip install -r requirements.txt
python analyze.py
python stats_analysis.py
```

No pinned versions are specified in `requirements.txt`; it just lists the packages imported by the scripts.

## Outputs

| Output file | Produced by | Description | Figure panel |
|---|---|---|---|
| `tmd-nanos_severity.pdf` | `analyze.py` | Stripplot + barplot (mean ± SD) of % severe-or-dead embryos (`dead` + `sev_malformed`, as a percentage of total scored) per condition, one point per repeat (n=3 repeats/condition). Y-axis: "% Severe or Dead"; x-axis: condition. | TODO |
| `statistical_results.txt` | `stats_analysis.py` | Text report: per-condition summary (n, mean, median, individual repeat values), a Kruskal-Wallis omnibus test across all four conditions, and pairwise Mann-Whitney U tests (U statistic, p-value, significance stars) between every pair of conditions. | TODO |

## Files

- `tmd-nanos_conditions.csv` — input data (see above).
- `analyze.py` — computes % severe-or-dead per condition/repeat and produces `tmd-nanos_severity.pdf`.
- `stats_analysis.py` — computes the same per-condition values and runs the statistical comparisons, writing `statistical_results.txt`.
- `utils.py` — shared helpers: `get_stars()` (p-value → significance stars) is used by `stats_analysis.py`. `parametric()` and `statistical_analysis()` are also defined here but are not currently called — with n=3 repeats per condition, `scipy.stats.normaltest` (which needs at least 8 samples) can't be used, so `stats_analysis.py` uses the non-parametric Mann-Whitney U / Kruskal-Wallis tests unconditionally instead.
