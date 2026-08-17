import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


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

print(severity_df)
# raise ValueError("Stop here for now")

color_palette = {
    "n_d": "#888888",
    "tmd-n_d": "#f5bb33",
    "wt": "#4a4a4a",
    "tmd-n_stuffer": "#c4c4c4",
}

sns.set_style("ticks")
plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 6,  # default text
        "axes.titlesize": 8,  # subplot titles
        "axes.labelsize": 6,  # x/y axis labels
        "axes.labelweight": "bold",  # x/y axis labels
        "xtick.labelsize": 6,  # x tick labels
        "ytick.labelsize": 6,  # y tick labels
        "legend.fontsize": 6,  # legend
        "figure.titlesize": 8,  # suptitle
    }
)

fig, ax = plt.subplots(figsize=(2, 2))

sns.stripplot(
    data=severity_df,
    y="pct_severe_or_dead",
    x="condition",
    hue="condition",
    # hue_order=["full_mix", "no_tdrd7a"],
    # palette=color_palette,
    palette="dark:black",
    # dodge=True,
    alpha=0.8,
    size=2,
    ax=ax,
    jitter=0.2,
)
# sns.pointplot(
#     data=severity_df,
#     y="pct_severe_or_dead",
#     hue="condition",
#     # hue_order=["full_mix", "no_tdrd7a"],
#     dodge=0.6,
#     ax=ax,
#     errorbar="sd",  # standard error
#     estimator="mean",  # or "mean"
#     capsize=0.075,
#     linestyle="none",
#     markersize=10,
#     marker="_",
#     err_kws=dict(linewidth=0.5, color="black"),
#     markeredgewidth=1,
#     palette="dark:black",
#     # zorder=5,
# )

sns.barplot(
    data=severity_df,
    y="pct_severe_or_dead",
    x="condition",
    hue="condition",
    # hue_order=["full_mix", "no_tdrd7a"],
    # dodge=True,
    ax=ax,
    errorbar="sd",  # standard error
    estimator="mean",  # or "mean"
    capsize=0.075,
    alpha=1,
    linewidth=0.5,
    err_kws=dict(linewidth=0.5, color="black"),
    edgecolor="black",
    palette=color_palette,
)

ax.set_ylabel("% Severe or Dead")
ax.set_xlabel("")
ax.set_title("TMD-Nanos Severity", fontweight="bold")
ax.set_box_aspect(1)


sns.despine()
plt.tight_layout()
fig.savefig("tmd-nanos_severity.pdf", transparent=True)
