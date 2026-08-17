from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from utils import statistical_analysis, get_stars


df = pd.read_csv(Path(__file__).parent / "data" / "nucleus.csv")

# Normalize the dnd and nanos fluorophore levels to the globin
df["dnd_normalized"] = df["dnd"] / df["globin"]
df["nanos_normalized"] = df["nanos"] / df["globin"]

# Normalize to full mix at each stage
df_stage_means = df.groupby(by=["stage", "condition"]).mean()
df["dnd_normalized_to_full_mix"] = df.apply(
    lambda row: (
        row["dnd_normalized"]
        / df_stage_means.loc[(row["stage"], "full_mix"), "dnd_normalized"]
    ),
    axis=1,
)
df["nanos_normalized_to_full_mix"] = df.apply(
    lambda row: (
        row["nanos_normalized"]
        / df_stage_means.loc[(row["stage"], "full_mix"), "nanos_normalized"]
    ),
    axis=1,
)


def format_p_value(p_value):
    if p_value < 0.0001:
        mantissa, exponent = f"{p_value:.2e}".split("e")
        return f"{mantissa}e{int(exponent)}"
    return f"{p_value:.4f}"


def plot_and_save_figures(df, swarm_plots=False):
    stats_lines = []

    # plt.rcParams["text.usetex"] = True
    sns.set_style("ticks")
    mpl.rcParams["mathtext.fontset"] = "stix"
    mpl.rcParams["font.family"] = "Arial"
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
    x_figure_size = 3.27  # 1 column figure size in inches for Cell Press journals
    y_figure_size = 1.63  # 1 column figure size in inches for Cell Press journals
    ncols = 2
    nrows = 2
    figure_name_list = ["normalized_graphs"]
    fig, ax = plt.subplots(
        nrows=nrows, ncols=ncols, figsize=(ncols * x_figure_size, nrows * y_figure_size)
    )

    err_kws_dict = {
        "capsize": 3,
        "elinewidth": 1,
        "marker": "D",
        "markersize": 4,
        "mew": 1,
        "mec": "#69696900",
        "barsabove": False,
    }
    line_style = ":"
    error_bar_param = "sd"

    # color_style_dict = {"full_mix": "#5a829e", "no_tdrd7": "#bfa35c"}
    color_style_dict = {
        "full_mix": "#888888",
        "no_tdrd7": "#7B3294",
    }
    if swarm_plots:
        figure_name_list.append("swarm_plots")
        for ax_idx, y_col in zip(
            [
                "dnd_normalized",
                "nanos_normalized",
                "dnd_normalized_to_full_mix",
                "nanos_normalized_to_full_mix",
            ],
            ax.flatten(),
        ):
            sns.stripplot(
                data=df,
                x="stage",
                y=ax_idx,
                ax=y_col,
                hue="condition",
                palette=color_style_dict,
                marker="o",
                jitter=0.12,
                linewidth=0.3,
                edgecolor="#050505",
                alpha=0.3,
                size=2.25,
            )

            stats_lines.append(f"Statistical analysis for {ax_idx}:")
            for stage in df["stage"].unique():
                stats_lines.append(f"\t {stage} stage:")
                data_set1 = (
                    df[(df["stage"] == stage) & (df["condition"] == "no_tdrd7")][ax_idx]
                    .dropna()
                    .to_list()
                )
                data_set2 = (
                    df[(df["stage"] == stage) & (df["condition"] == "full_mix")][ax_idx]
                    .dropna()
                    .to_list()
                )
                test_type, _, p_value = statistical_analysis(data_set1, data_set2)
                stats_lines.append(
                    f"\t \t test-type = {test_type}{' ' * (19 - (int(len(test_type))))}, p-value = {format_p_value(p_value)}, significance = {get_stars(p_value)}, n(no_tdrd7) = {len(data_set1)}, n(full_mix) = {len(data_set2)}"
                )

    sns.lineplot(
        data=df,
        x="stage",
        y="dnd_normalized",
        ax=ax[0, 0],
        hue="condition",
        err_style="bars",
        errorbar=error_bar_param,
        err_kws=err_kws_dict,
        linestyle=line_style,
        palette=color_style_dict,
    )

    sns.lineplot(
        data=df,
        x="stage",
        y="nanos_normalized",
        ax=ax[0, 1],
        hue="condition",
        err_style="bars",
        errorbar=error_bar_param,
        err_kws=err_kws_dict,
        linestyle=line_style,
        palette=color_style_dict,
    )

    sns.lineplot(
        data=df,
        x="stage",
        y="dnd_normalized_to_full_mix",
        ax=ax[1, 0],
        hue="condition",
        err_style="bars",
        errorbar=error_bar_param,
        err_kws=err_kws_dict,
        linestyle=line_style,
        palette=color_style_dict,
    )
    ax[1, 0].set_ylim(None, 1.6)
    ax[1, 0].set_xlabel("Stage")
    ax[1, 0].set_title("Normailzed mean nuclear fluorescence across stages")

    sns.lineplot(
        data=df,
        x="stage",
        y="nanos_normalized_to_full_mix",
        ax=ax[1, 1],
        hue="condition",
        err_style="bars",
        errorbar=error_bar_param,
        err_kws=err_kws_dict,
        linestyle=line_style,
        palette=color_style_dict,
    )
    ax[1, 1].set_ylim(None, 2.4)
    ax[1, 1].set_xlabel("Stage")
    ax[1, 1].set_title("Normailzed nuclear fluorescence across stages")

    sns.despine()
    plt.tight_layout()
    figure_name_list.append(".pdf")
    fig_name = "_".join(figure_name_list)
    with mpl.rc_context({"mathtext.fontset": "stix", "font.family": "STIXGeneral"}):
        ax[1, 0].set_ylabel(
            r"$\frac{\overline{I}_{eGFP.dnd1 3'UTR}}{\overline{I}_{tagBFP.globin 3'UTR}}$"
        )
        ax[1, 1].set_ylabel(
            r"$\frac{\overline{I}_{mScarlet-i.nos3 3'UTR}}{\overline{I}_{tagBFP.globin 3'UTR}}$"
        )

    fig.savefig(Path(__file__).parent / "output" / fig_name)

    with open(Path(__file__).parent / "output" / "stats.txt", "w") as f:
        print(len(stats_lines))
        for line in stats_lines:
            f.write(line + "\n")


plot_and_save_figures(df, swarm_plots=True)
# plot_and_save_figures(df, swarm_plots=False)
