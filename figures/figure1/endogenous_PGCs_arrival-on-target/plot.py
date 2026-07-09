from numpy import test
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from pathlib import Path
from utils import statistical_analysis, get_stars

root = Path(__file__).parent
fig_dir = root / "figure"
fig_dir.mkdir(exist_ok=True)

df = pd.read_csv(root / "Tdrd7a_ectopicity.csv")

df["total_number"] = df["Gonad"] + df["Ectopic"]
df["ectopicity"] = df["Ectopic"] / df["total_number"] * 100

nrows = 1
ncols = 2

width_ax = 1.1
height_ax = 1.8

sns.set_style("ticks")

plt.rcParams.update(
    {
        "font.size": 6,  # default text
        "font.family": "Arial",
        "axes.titlesize": 8,  # subplot titles
        "axes.labelsize": 6,  # x/y axis labels
        "xtick.labelsize": 6,  # x tick labels
        "ytick.labelsize": 6,  # y tick labels
        "legend.fontsize": 6,  # legend
        "figure.titlesize": 8,  # suptitle
    }
)

fig, ax = plt.subplots(
    ncols=ncols, nrows=nrows, figsize=(ncols * width_ax, nrows * height_ax)
)
color_palette = {
    "ctrl": "#878787",
    "kd": "#7B3294",
}

stats_list = []
for ax_idx, feature in enumerate(["ectopicity", "total_number"]):
    # sns.stripplot(
    #    data=df,
    #    y=feature,
    #    hue="condition",
    #    palette=color_palette,
    #    hue_order=list(color_palette.keys()),
    #    dodge=0.4,
    #    ax=ax[ax_idx],
    #    alpha=0.4,
    #    size=2.5,
    # )

    # sns.stripplot(
    #    data=df,
    #   y=feature,
    #    hue="condition",
    #    palette="dark:black",
    #    alpha=0.7,
    #    ax=ax[ax_idx],
    #    dodge=0.4,
    #    size=2,
    # )

    sns.stripplot(
        data=df,
        y=feature,
        hue="condition",
        hue_order=list(color_palette.keys())
        if isinstance(color_palette, dict)
        else None,
        palette=color_palette,
        dodge=0.4,
        alpha=0.4,
        size=2,
        ax=ax[ax_idx],
        jitter=0.3,
    )

    sns.pointplot(
        data=df,
        y=feature,
        hue="condition",
        hue_order=list(color_palette.keys())
        if isinstance(color_palette, dict)
        else None,
        dodge=0.4,
        ax=ax[ax_idx],
        errorbar="sd",  # standard error
        estimator="median",  # or "mean"
        capsize=0.075,
        linestyle="none",
        markersize=10,
        marker="_",
        err_kws=dict(linewidth=0.5, color="black"),
        markeredgewidth=1,
        palette="dark:black",
        zorder=5,
    )

    if feature == "ectopicity":
        ax[ax_idx].set_ylim(-1, None)
    ax[ax_idx].set_title(feature, fontweight="bold")

    # Statistical Analysis
    stats_list.append(f"{feature}:\n")
    data_ctrl = df[df["condition"] == "ctrl"][feature].to_list()
    data_kd = df[df["condition"] == "kd"][feature].to_list()

    test_type, _, p_value = statistical_analysis(data_ctrl, data_kd)
    stats_list.extend(
        [
            f"\t Test:    {test_type}\n"
            f"\t p-Value: {p_value}\n"
            f"\t stars:   {get_stars(p_value)}\n",
            f"\t n ctrl:  {len(data_ctrl)}\n",
            f"\t n KD:  {len(data_ctrl)}\n",
        ]
    )


with open(fig_dir / "stats.txt", "w") as f:
    f.writelines(stats_list)

sns.despine()
plt.tight_layout()
fig.savefig(fig_dir / "ePGCs_ectopicity_number.pdf")
plt.show()
