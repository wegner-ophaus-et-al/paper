import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
import seaborn as sns
from pathlib import Path
from scipy.stats import gaussian_kde
import numpy as np

mpl.rcParams["font.family"] = "Arial"


def contact_sheet(list_of_cells: list, save_dir: Path):
    """
    Generate a contanct sheet for all cells in the list
    """
    fig, axs = plt.subplots(len(list_of_cells), 3, figsize=(10, 4 * len(list_of_cells)))
    for ax, cell in zip(axs, list_of_cells):
        cell.plot_markers_on_axis(ax, segmentation_cmap="summer")
        for a in ax:
            a.axis("off")

    figure_dir = save_dir / "figures"
    figure_dir.mkdir(exist_ok=True)
    fig.savefig(figure_dir / "contact_sheet.pdf")


def plot_individial_granule_profile(
    df, save_dir: Path, color_palette: str | dict = "Set2"
):
    """
    Plot the granule profile for all cells in the dataframe
    """
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
    fig = plt.figure(figsize=(12, 6))
    fig.suptitle("Granule Profile for Individual Granules", fontsize=8)

    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1])

    axs = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
    ]

    for ax, feature in zip(
        axs,
        [
            "Area",
            "EdgeDistanceToNucleus",
            "CentroidDistanceToNucleus",
            "GranuleKurtosis",
        ],
    ):
        sns.violinplot(
            data=df,
            x="marker",
            y=feature,
            hue="condition",
            palette=color_palette,
            hue_order=list(color_palette.keys())
            if isinstance(color_palette, dict)
            else None,
            ax=ax,
            split=True,
            inner="quart",
            legend=False,
        )
    fig.tight_layout()
    sns.despine()
    figure_dir = save_dir / "figures"
    figure_dir.mkdir(exist_ok=True)
    fig.savefig(figure_dir / "Individual.pdf")


def per_cell_summary(
    df, save_dir: Path, aggregation_function="mean", color_palette: str | dict = "Set2"
):
    """
    Aggregate the granule features per cell and plot the summary statistics
    """
    df_agg = (
        df.groupby(["uid", "marker", "condition", "stage"])
        .agg(aggregation_function, numeric_only=True)
        .reset_index()
    )

    df_agg_sum = (
        df.groupby(["uid", "marker", "condition", "stage"])
        .agg("sum", numeric_only=True)
        .reset_index()
    )

    df_agg_std = (
        df.groupby(["uid", "marker", "condition", "stage"])
        .agg("std", numeric_only=True)
        .reset_index()
    )

    features_mean = [
        "Area",
        "SphericalVolume",
        "GranuleNumberPerCell",
        "EdgeDistanceToNucleusSigned",
        "CentroidDistanceToNucleusSigned",
        "TouchAreaNucleus",
    ]

    features_sum = ["Area", "SphericalVolume", "EllipsoidVolumeProlate"]

    features_std = ["Area", "SphericalVolume", "GranuleNumberPerCell"]

    ncols = 6
    nrows = 4
    subfig_cw = 1.2
    subfig_ch = 2

    fig = plt.figure(figsize=(ncols * subfig_cw, subfig_ch * nrows))

    gs = gridspec.GridSpec(ncols=ncols, nrows=nrows)

    axs = [
        [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 3])],
        [fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 4])],
        [fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[0, 5])],
        [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 3])],
        [fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 4])],
        [fig.add_subplot(gs[1, 2]), fig.add_subplot(gs[1, 5])],
    ]
    axs_sum = [
        [fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 3])],
        [fig.add_subplot(gs[2, 1]), fig.add_subplot(gs[2, 4])],
        [fig.add_subplot(gs[2, 2]), fig.add_subplot(gs[2, 5])],
    ]

    axs_std = [
        [fig.add_subplot(gs[3, 0]), fig.add_subplot(gs[3, 3])],
        [fig.add_subplot(gs[3, 1]), fig.add_subplot(gs[3, 4])],
        [fig.add_subplot(gs[3, 2]), fig.add_subplot(gs[3, 5])],
    ]

    log_scale_features = [
        "Area",
        "SphericalVolume",
        "EllipsoidVolumeProlate",
        # "EdgeDistanceToNucleusSigned",
    ]

    clip_to_ten_features = [
        "EdgeDistanceToNucleusSigned",
        "CentroidDistanceToNucleusSigned",
    ]

    for feature, ax in zip(features_mean, axs):
        for ax_idx, marker in enumerate(["gra", "dnd1"]):
            if feature in log_scale_features:
                ax[ax_idx].set_yscale("log")
            if feature in clip_to_ten_features:
                df_agg[feature] = df_agg[feature].clip(-10, 8)
            sns.violinplot(
                data=df_agg[df_agg["marker"] == marker],
                x="stage",
                y=feature,
                order=["8hpf", "24hpf"],
                hue="condition",
                hue_order=list(color_palette.keys())
                if isinstance(color_palette, dict)
                else None,
                palette=color_palette,
                split=True,
                inner="quart",
                ax=ax[ax_idx],
                linewidth=0.8,
                legend=False,
            )
            ax[ax_idx].set_title(f"Title", fontsize=8, fontweight="bold")

    for feature, ax in zip(features_sum, axs_sum):
        for ax_idx, marker in enumerate(["gra", "dnd1"]):
            if feature in log_scale_features:
                ax[ax_idx].set_yscale("log")
            sns.violinplot(
                data=df_agg_sum[df_agg_sum["marker"] == marker],
                x="stage",
                y=feature,
                hue="condition",
                hue_order=list(color_palette.keys())
                if isinstance(color_palette, dict)
                else None,
                palette=color_palette,
                split=True,
                inner="quart",
                ax=ax[ax_idx],
                linewidth=0.8,
                legend=False,
            )
            ax[ax_idx].set_title(f"Title", fontsize=8, fontweight="bold")

    for feature, ax in zip(features_std, axs_std):
        for ax_idx, marker in enumerate(["gra", "dnd1"]):
            if feature in log_scale_features:
                ax[ax_idx].set_yscale("log")
            sns.violinplot(
                data=df_agg_std[df_agg_std["marker"] == marker],
                x="stage",
                y=feature,
                hue="condition",
                hue_order=list(color_palette.keys())
                if isinstance(color_palette, dict)
                else None,
                palette=color_palette,
                split=True,
                inner="quart",
                ax=ax[ax_idx],
                linewidth=0.8,
                legend=False,
            )
            ax[ax_idx].set_title(f"Title", fontsize=8, fontweight="bold")

    plt.tight_layout()
    sns.despine()
    figure_dir = save_dir / "figures"
    figure_dir.mkdir(exist_ok=True)
    fig.savefig(
        figure_dir / f"plot_per_granule_profile_aggregate-{aggregation_function}.pdf",
        transparent=True,
    )


def plot_foldchange(df, save_dir: Path, color_palette="Set2"):
    """
    Plot the differences between conditions for each marker and stage
    """
    # Aggregate per marker, condition and cell by mean
    df_agg = (
        df.groupby(["marker", "condition", "stage"])
        .agg("mean", numeric_only=True)
        .reset_index()
    )

    stage_dfs = {
        "8hpf": df_agg[df_agg["stage"] == "8hpf"],
        "24hpf": df_agg[df_agg["stage"] == "24hpf"],
    }

    # Get fold change between conditions for each marker per stage
    for stage, stage_df in stage_dfs.items():
        fig, ax = plt.subplots(figsize=(6, 4))
        for marker in stage_df["marker"].unique():
            marker_df = stage_df[stage_df["marker"] == marker]
            fold_change = (
                marker_df[marker_df["condition"] == "kd"]["Area"].values[0]
                / marker_df[marker_df["condition"] == "ctrl"]["Area"].values[0]
            )
            ax.bar(marker, fold_change, color=color_palette[0], alpha=0.7)

        ax.set_title(f"Fold Change (KD / CTRL) for {stage}", fontsize=10)
        ax.set_ylabel("Fold Change (Area)", fontsize=8)
        sns.despine()
        figure_dir = save_dir / "figures"
        figure_dir.mkdir(exist_ok=True)
        fig.savefig(figure_dir / f"fold_change_{stage}.pdf")


def plot_ridgeplot(df, marker: str, ax, feature: str = "featureX"):
    """
    Draw a ridge plot for a given marker onto an existing Axes.

    Groups by (stage, condition) — each combination gets its own KDE ridge.

    Parameters
    ----------
    df      : DataFrame with columns [stage, condition, marker, featureX]
    marker  : which marker value to filter on
    ax      : the matplotlib Axes to draw into
    feature : column name to plot (default "featureX")
    """
    sub = df[df["marker"] == marker].copy()

    # Build one ridge per (stage, condition) group
    groups = sub.groupby(["stage", "condition"])
    group_keys = list(groups.groups.keys())
    n = len(group_keys)

    # Color palette — one color per group
    colors = sns.color_palette("Set2", n_colors=n)

    # Global x range for consistent KDE evaluation
    xmin, xmax = sub[feature].min(), sub[feature].max()
    xs = np.linspace(xmin, xmax, 300)

    # Overlap factor: fraction of vertical space each ridge occupies
    overlap = 1.4
    row_height = 1.0 / n

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(0, 1)
    ax.axis("off")  # hide the outer frame; we draw everything via inset axes

    for i, (key, color) in enumerate(zip(group_keys, colors)):
        stage, condition = key
        vals = groups.get_group(key)[feature].dropna()

        # KDE
        kde = gaussian_kde(vals, bw_method=0.3)
        ys = kde(xs)
        ys = ys / ys.max()  # normalize to [0, 1]

        # Position this ridge: bottom edge in axis-fraction coordinates
        bottom = 1.0 - (i + 1) * row_height
        height = row_height * overlap

        inset = ax.inset_axes([0, bottom, 1, height])  # [x0, y0, w, h] in ax fractions
        inset.set_xlim(xmin, xmax)
        inset.set_ylim(-0.05, 1.1)
        inset.patch.set_alpha(0)  # transparent background so ridges overlap nicely

        # Filled KDE
        inset.fill_between(xs, ys, alpha=0.8, color=color)
        # White outline on top
        inset.plot(xs, ys, color="white", lw=1.5)
        # Baseline
        inset.axhline(0, color=color, lw=1.5, clip_on=False)

        inset.axis("off")

        # Label on the left
        inset.text(
            -0.01,
            0.1,
            f"{stage}\n{condition}",
            transform=inset.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            fontweight="bold",
            color=color,
        )

    ax.set_title(f"Marker: {marker}", fontsize=11, fontweight="bold", pad=4)


def ridgeplot_per_marker(df, save_dir: Path, color_palette="Set2"):
    """
    Plot ridgeplots for each marker and stage
    """
    markers = df["marker"].unique()
    stages = df["stage"].unique()

    fig, axs = plt.subplots(len(stages), len(markers), figsize=(12, 6))

    for i, stage in enumerate(stages):
        for j, marker in enumerate(markers):
            plot_ridgeplot(df[df["stage"] == stage], marker, axs[i, j], feature="Area")

    plt.tight_layout()
    sns.despine()
    figure_dir = save_dir / "figures"
    figure_dir.mkdir(exist_ok=True)
    fig.savefig(figure_dir / "ridgeplots_per_marker.pdf")
