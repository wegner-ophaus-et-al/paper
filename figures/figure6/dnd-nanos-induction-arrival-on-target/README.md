# Quantification of Dnd1 + Nanos3 induction of germ cell fate

Zebrafish embryos were injected with either a full germ-cell-fate induction mRNA mix (dnd1, nanos3, buc, tdrd7a, tdrd6, vasa — condition `Full mix`) or a minimal mix of just dnd1 + nanos3 (condition `DndNos`). At 24 hpf, induced primordial germ cells (iPGCs) were scored per embryo as localizing to the gonad, showing high ligand signal, or being ectopic (mislocalized). This analysis computes the fraction of iPGCs "on target" (gonad + high-ligand) per embryo and compares the two conditions statistically.

## Input data

- `iPGC_localization_24hpf_FullmixvsDndNos.csv` — one row per embryo, columns `Condition` (`Full mix` / `DndNos`), `Gonad`, `LigandHigh`, `Ectopic` (cell counts), `Date` (replicate date). 80 rows total (2 conditions × 2 replicate dates × 20 embryos).
- The CSV was exported/fused from a source `.xlsx` workbook (per-replicate sheets plus the combined sheet used for the CSV); the workbook itself is not tracked in this repo (`*.xlsx` is gitignored). The CSV is the rawest form of the data included here; no upstream imaging data is included in this repo.

## Outputs

- `dnd-nos_induced-iPGCs_ectopicity.pdf` — strip plot of on-target iPGC percentage per embryo by condition, with median ± SD overlay.
- `iPGC_dnd_nanos_ectopicity_stats.txt` — test used (t-test or Mann-Whitney U, chosen via a normality check), medians, means, counts, and p-value/significance for the `Full mix` vs `DndNos` comparison.

## Key files

- `analyze_and_plot.py` — loads the CSV, computes the on-target percentage, generates the plot, and runs/writes the statistics.
- `utils.py` — statistics helpers: normality check, test selection (t-test vs. Mann-Whitney U), and significance-star formatting.
