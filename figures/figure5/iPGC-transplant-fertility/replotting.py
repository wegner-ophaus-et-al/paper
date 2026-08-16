import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
from pathlib import Path


sns.set_style("ticks")

mpl.rcParams["font.family"] = "Arial"

plt.rcParams.update(
    {
        "font.size": 6,  # default text
        "font.family": "Arial",
        "axes.titlesize": 8,  # subplot titles
        "axes.labelsize": 6,  # x/y axis labels
        "axes.labelweight": "bold",  # x/y axis labels
        "xtick.labelsize": 6,  # x tick labels
        "ytick.labelsize": 6,  # y tick labels
        "legend.fontsize": 6,  # legend
        "figure.titlesize": 8,  # suptitle
    }
)

root = Path(__file__).parent
df = pd.read_csv(root / "data" / "iPGC_transplantation_fertility_data.csv")

df["percent_fertilized"] = (
    df["fertilized"] / (df["fertilized"] + df["not_fertilized"])
) * 100

df["percent_dominant_marker"] = (
    df["red_eye_positive"] / (df["red_eye_positive"] + df["red_eye_negative"])
) * 100

# Kick out the fish that have dominant maker <100%, but keep the empty/nan
df = df[df["percent_dominant_marker"].isna() | (df["percent_dominant_marker"] == 100)]

df["is_fertile"] = df["percent_fertilized"] > 0
df["is_dominant_marker"] = df["percent_dominant_marker"] > 0

df.to_csv(
    root / "data" / "iPGC_transplantation_fertility_data_processed.csv", index=False
)


df.groupby("condition").describe(include="all").to_csv(
    root / "data" / "iPGC_transplantation_fertility_data_summary.csv"
)

color_palette = {
    "wt": "#444444",
    "full_mix": "#888888",
    "no_tdrd7a": "#7B3294",
}


ncols = 5
nrows = 1
width = 1.4
height = 2.5

fig, axs = plt.subplots(
    nrows=nrows, ncols=ncols, figsize=(ncols * width, nrows * height)
)

for idax, feature in enumerate(
    [
        "percent_fertilized",
        "percent_dominant_marker",
    ]
):
    sns.stripplot(
        data=df[df["condition"] != "wt"],
        y=feature,
        hue="condition",
        hue_order=["full_mix", "no_tdrd7a"],
        palette=color_palette,
        dodge=0.4,
        alpha=0.4,
        size=2,
        ax=axs[idax],
        jitter=0.2,
    )
    sns.pointplot(
        data=df[df["condition"] != "wt"],
        y=feature,
        hue="condition",
        hue_order=["full_mix", "no_tdrd7a"],
        dodge=0.4,
        ax=axs[idax],
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
    title_map = {
        "percent_fertilized": "Fertilization rate (%)",
        "percent_dominant_marker": "Germline contribution (%)",
    }
    axs[idax].set_title(title_map[feature])

    axs[idax].get_legend().remove()  # remove legend for bar plots
for axid, feature in enumerate(["is_fertile", "is_dominant_marker"]):
    sns.barplot(
        data=df if feature == "is_fertile" else df[df["condition"] != "wt"],
        x="condition",
        y=feature,
        hue="condition",
        order=["wt", "full_mix", "no_tdrd7a"]
        if feature == "is_fertile"
        else ["full_mix", "no_tdrd7a"],
        errorbar=("ci", 95),
        ax=axs[axid + 2],
        capsize=0.075,
        palette=color_palette,
        err_kws=dict(linewidth=0.5, color="black"),
    )
    title_map = {
        "is_fertile": "Proportion fertile",
        "is_dominant_marker": "Proportion with dominant marker",
    }
    axs[axid + 2].set_title(title_map[feature])
    if feature == "is_fertile":
        axs[axid + 2].yaxis.set_major_formatter(
            plt.FuncFormatter(lambda y, _: f"{y * 100:.0f}")
        )
        axs[axid + 2].set_ylabel("Proportion fertile (%)")


# Sex ratios
sex_ratios = {
    "full_mix": float(
        len(df[df["condition"] == "full_mix"]) / len(df[df["condition"] == "full_mix"])
    ),  # Since all fish developed males count
    "no_tdrd7a": float(
        len(df[df["condition"] == "no_tdrd7a"])
        / len(df[df["condition"] == "no_tdrd7a"])
    ),  # Since all fish developed males count
    "wt": float(24 / (24 + 18)),  # 24 males and 18 females
}

sns.barplot(
    x=sex_ratios.keys(),
    y=sex_ratios.values(),
    hue=sex_ratios.keys(),
    order=["wt", "full_mix", "no_tdrd7a"],
    palette=color_palette,
    ax=axs[4],
)
axs[4].set_title("Sex ratio")
axs[4].yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y * 100:.0f}"))
axs[4].set_ylabel("Percent of males on clutch [%]")


sns.despine()
plt.tight_layout()
fig_save_path = root / "figures" / "iPGC_transplantation_fertility_replot.pdf"
fig_save_path.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(fig_save_path)
plt.show()
