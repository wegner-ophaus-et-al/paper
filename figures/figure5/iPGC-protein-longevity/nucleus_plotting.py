from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl


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


def plot_and_save_figures(df, swarm_plots=False):
    # plt.rcParams["text.usetex"] = True
    mpl.rcParams["mathtext.fontset"] = "stix"
    mpl.rcParams["font.family"] = "STIXGeneral"
    sns.set_style("ticks")
    figure_name_list = ["normalized_graphs"]
    fig, ax = plt.subplots(2, 2, figsize=(8, 5))

    err_kws_dict = {
        "capsize": 5,
        "marker": "D",
        "mew": 1,
        "mec": "white",
        "barsabove": True,
    }
    line_style = ":"
    error_bar_param = "sd"

    color_style_dict = {"full_mix": "#5a829e", "no_tdrd7": "#bfa35c"}

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
                # dodge=True,
                # palette="dark:.25",
                palette=color_style_dict,
                alpha=0.5,
                size=2,
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


plot_and_save_figures(df, swarm_plots=True)
plot_and_save_figures(df, swarm_plots=False)
