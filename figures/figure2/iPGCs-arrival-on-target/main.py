import marimo

__generated_with = "0.23.9"
app = marimo.App()


@app.cell
def _():
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    from utils import statistical_analysis, get_stars

    return get_stars, pd, plt, sns, statistical_analysis


@app.cell
def _(pd):
    df_raw = pd.read_csv("ectopicity_data.csv")
    df_raw["EmbryoID"] = df_raw.index
    return (df_raw,)


@app.cell
def _(df_raw, pd):
    df = pd.melt(
        df_raw,
        id_vars=["Date", "Condition", "EmbryoID"],  # Keep only the identifying columns
        value_vars=["Gonad", "Ectopic", "HighLigandArea"],
        var_name="Region",
        value_name="Count",
    )
    df["Region_clean"] = df["Region"].copy()
    df.loc[df["Region_clean"] == "HighLigandArea", "Region_clean"] = "Gonad"
    df.loc[df["Region_clean"] == "OnTarget", "Region_clean"] = "Gonad"
    return (df,)


@app.cell
def _():
    MIN_CELL_NUMBER = 10
    palette_dict = {
        "with Tdrd7a": "#1d8d91",
        "without Tdrd7a": "#E99949",
    }
    return (palette_dict,)


@app.cell
def _(df):
    # Divide the Count by the total number of counts per EmbryoID
    df["TotalCount"] = df.groupby(["Date", "Condition", "EmbryoID"])["Count"].transform(
        "sum"
    )
    df["Percent"] = (df["Count"] / df["TotalCount"]) * 100
    return


@app.cell
def _(df):
    df
    return


@app.cell
def _(df, palette_dict, plt, sns):
    # Set font to SF Pro Semibold
    plt.rcParams["font.family"] = "Arial"
    # increase the font size and weight

    # sns.set_context("poster", font_scale=1.5)

    fig_size = 10
    fig_aspect = 1.4
    sns.set_style("ticks")
    fig, ax = plt.subplots(figsize=(fig_size, fig_aspect * fig_size))
    sns.violinplot(
        data=df.query("TotalCount > @MIN_CELL_NUMBER"),
        y="Percent",
        x="Region",
        hue="Condition",
        split=True,
        inner="quart",
        palette=palette_dict,
        ax=ax,
        linewidth=4,
        dodge=True,
        bw_method=0.4,
    )
    ax.set_ylabel("Percent of cell in region")
    ax.set_xlabel("Region")

    sns.despine()
    fig.savefig("ectopicity_violinplot.pdf", transparent=True, bbox_inches="tight")
    return


@app.cell
def _(df):
    df.query("Region == 'Ectopic'")
    return


@app.cell
def _(df, plt, sns):
    sns.set_style("ticks")
    plt.rcParams.update(
        {
            "font.size": 6,  # default text
            "axes.titlesize": 8,  # subplot titles
            "axes.labelsize": 6,  # x/y axis labels
            "xtick.labelsize": 6,  # x tick labels
            "ytick.labelsize": 6,  # y tick labels
            "legend.fontsize": 6,  # legend
            "figure.titlesize": 8,  # suptitle
        }
    )
    plt.rcParams["font.family"] = "Arial"

    color_palette = {
        "with Tdrd7a": "#888888",
        "without Tdrd7a": "#7B3294",
    }
    # sns.set_context(context="paper")

    fig_paper, ax_paper = plt.subplots(1, 1, figsize=(1.35, 2))

    sns.stripplot(
        data=df.query("Region == 'Ectopic'"),
        x="Region",
        y="Percent",
        hue="Condition",
        # hue_order=list(color_palette.keys()),
        palette=color_palette,
        dodge=0.4,
        alpha=0.4,
        size=3,
        ax=ax_paper,
        jitter=0.2,
    )

    sns.pointplot(
        data=df.query("Region == 'Ectopic'"),
        x="Region",
        y="Percent",
        hue="Condition",
        dodge=0.4,
        ax=ax_paper,
        errorbar="sd",  # or 'se' standard error
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
    ax_paper.set_xlabel("")
    ax_paper.set_xticks([])
    ax_paper.set_xticklabels([])
    ax_paper.set_ylabel("Ectopic iPGCs [%]", fontweight="bold")
    ax_paper.set_title("Ectopic iPGCs at 24 hpf", fontweight="bold")
    ax_paper.legend_.remove()
    plt.tight_layout()
    sns.despine()
    plt.show()
    fig_paper.savefig("ectopicity_paper_plot.pdf")
    return


@app.cell
def _(df):
    for condition in df["Condition"].unique():
        count = len(
            df.query("TotalCount > @MIN_CELL_NUMBER")
            .groupby(["Condition", "EmbryoID"])
            .count()
            .query("Condition == @condition")
        )
        print(f"For condition <{condition}> the sample size is {count}")
    return


@app.cell
def _(df):
    df.groupby(["Condition", "Region"])["Percent"].median()
    return


@app.cell
def _(df):
    df
    return


@app.cell
def _(df, get_stars, statistical_analysis):
    ctrl_values = df.loc[
        (df["Condition"] == "with Tdrd7a") & (df["Region"] == "Ectopic"), "Percent"
    ].to_list()
    nt_values = df.loc[
        (df["Condition"] == "without Tdrd7a") & (df["Region"] == "Ectopic"), "Percent"
    ].to_list()

    test, _, p_value = statistical_analysis(ctrl_values, nt_values)
    stars = get_stars(p_value)
    print(
        f"Ectopic PGC between full mix and mix without tdrd7a - {test}:\n",
        f"p value: {round(p_value, 6)} \t -> {stars}\n",
        f"Count FM: {len(ctrl_values)}\n",
        f"Count NT: {len(nt_values)}",
    )
    return ctrl_values, nt_values, p_value, stars, test


@app.cell
def _(ctrl_values, nt_values, p_value, stars, test):
    lines = [
        "Ectopic PGC between full mix and mix without tdrd7a\n",
        f"\t Count full mix: {len(ctrl_values)}\n",
        f"\t Count w/o tdrd: {len(nt_values)}\n",
        f"\t Test type:      {test}\n",
        f"\t p value:        {round(p_value, 5)}\n",
        f"\t stars:          {stars}\n",
    ]
    return (lines,)


@app.cell
def _(lines):
    with open("stats.txt", "w+") as f:
        f.writelines(lines)
    return


@app.cell
def _():

    return


if __name__ == "__main__":
    app.run()
