import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

df = pd.read_csv(
    "/Users/julian/Library/Mobile Documents/com~apple~CloudDocs/Documents/General Science/Programming/py/general_analysis/20260807_nanos_dnd/iPGC_localization_24hpf_FullmixvsDndNos.csv"
)
df["OnTarget"] = (
    (df["Gonad"] + df["LigandHigh"]) / (df["Ectopic"] + df["LigandHigh"] + df["Gonad"])
) * 100
print(df.head())

palette = {"Full mix": "#888888", "DndNos": "#FFCD22"}

sns.set_theme(style="ticks")
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
fig, ax = plt.subplots(figsize=(1.6, 2.3))

sns.stripplot(
    data=df,
    y="OnTarget",
    x="Condition",
    hue="Condition",
    alpha=0.4,
    size=3,
    ax=ax,
)
sns.pointplot(
    data=df,
    x="Condition",
    y="OnTarget",
    hue="Condition",
    # dodge=0.4,
    errorbar="sd",  # standard error
    estimator="median",  # or "median"
    capsize=0.075,
    linestyle="none",
    markersize=10,
    marker="_",
    err_kws=dict(linewidth=0.4, color="black"),
    markeredgewidth=1,
    palette="dark:black",
    zorder=5,
    ax=ax,
)

ax.set_ylim(50, 105)
ax.set_title("On target Dnd1 + Nanos3 induced PGCs")
sns.despine()
plt.tight_layout()
plt.show()
fig.savefig("dnd-nos_induced-iPGCs_ectopicity.pdf", transparent=True)
