from itertools import combinations

import pandas as pd
from scipy import stats

from utils import get_stars

df = pd.read_csv("tmd-nanos_conditions.csv")

grouped = df.groupby(["condition", "repeat"])
totals = grouped["count"].sum().rename("total")
bad = (
    df[df["category"].isin(["dead", "sev_malformed"])]
    .groupby(["condition", "repeat"])["count"]
    .sum()
    .rename("bad_count")
)

severity_df = pd.concat([totals, bad], axis=1).fillna(0)
severity_df["pct_severe_or_dead"] = (
    100 * severity_df["bad_count"] / severity_df["total"]
)
severity_df = severity_df.reset_index()

conditions = sorted(severity_df["condition"].unique())
values_by_condition = {
    cond: severity_df.loc[
        severity_df["condition"] == cond, "pct_severe_or_dead"
    ].tolist()
    for cond in conditions
}

lines = []
lines.append("Statistical analysis: % severe or dead by condition")
lines.append("=" * 55)
lines.append("")
lines.append(
    "Each condition has n=3 repeats. scipy.stats.normaltest (used by "
    "utils.parametric) requires at least 8 samples per group and raises "
    "an error below that, so normality cannot be formally assessed here. "
    "All comparisons below therefore use the non-parametric Mann-Whitney U "
    "test unconditionally, and an omnibus Kruskal-Wallis test is used in "
    "place of a one-way ANOVA."
)
lines.append(
    "Caveat: with n=3 vs n=3, the exact two-sided Mann-Whitney U test has "
    "a minimum possible p-value of 0.1, so no pairwise comparison here can "
    "reach the conventional p<0.05 threshold. Exact p-values are reported "
    "below rather than relying on significance stars alone."
)
lines.append("")

lines.append("Per-condition summary")
lines.append("-" * 55)
for cond in conditions:
    vals = values_by_condition[cond]
    n = len(vals)
    mean = sum(vals) / n
    median = stats.scoreatpercentile(vals, 50)
    lines.append(
        f"{cond:>15}: n={n}, mean={mean:6.2f}%, median={median:6.2f}%, "
        f"values={[round(v, 2) for v in vals]}"
    )
lines.append("")

lines.append("Omnibus test across all conditions")
lines.append("-" * 55)
h_statistic, kw_p = stats.kruskal(*[values_by_condition[c] for c in conditions])
lines.append(
    f"Kruskal-Wallis H-test: H={h_statistic:.4f}, p={kw_p:.4f} {get_stars(kw_p)}"
)
lines.append("")

lines.append("Pairwise comparisons (Mann-Whitney U test)")
lines.append("-" * 55)
lines.append(f"{'condition A':>15} {'condition B':>15} {'U':>8} {'p-value':>10}  stars")
for cond_a, cond_b in combinations(conditions, 2):
    u_statistic, p_value = stats.mannwhitneyu(
        values_by_condition[cond_a], values_by_condition[cond_b]
    )
    stars = get_stars(p_value)
    lines.append(
        f"{cond_a:>15} {cond_b:>15} {u_statistic:8.2f} {p_value:10.4f}  {stars}"
    )

output_text = "\n".join(lines) + "\n"

with open("statistical_results.txt", "w") as f:
    f.write(output_text)

print(output_text)
