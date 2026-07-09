import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
import matplotlib.lines as mlines
import seaborn as sns
from pathlib import Path
from scipy.stats import gaussian_kde
import numpy as np
from utils import statistical_analysis, get_stars

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
    plt.close(fig)


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
    ncols = 3
    nrows = 2
    w_size = 1.63
    h_size = 3.27
    fig = plt.figure(figsize=(ncols * w_size, nrows * h_size))
    fig.suptitle("Granule Profile for Individual Granules", fontsize=8)

    gs = gridspec.GridSpec(nrows=nrows, ncols=ncols)

    axs = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
        fig.add_subplot(gs[0, 2]),
        fig.add_subplot(gs[1, 2]),
    ]

    for ax, feature in zip(
        axs,
        [
            "Area",
            "EdgeDistanceToNucleus",
            "CentroidDistanceToNucleus",
            "GranuleKurtosis",
            "RelativeNuclearDistance",
            "SphericalVolume",
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
    figure_dir = save_dir / "figures" / "individial_granule_data"
    figure_dir.mkdir(exist_ok=True, parents=True)
    fig.savefig(figure_dir / "Individual.pdf")
    plt.close(fig)


def per_cell_summary(
    df, save_dir: Path, aggregation_function="mean", color_palette: str | dict = "Set2"
):
    """
    Aggregate the granule features per cell and plot the summary statistics
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
    nrows = 6
    subfig_cw = 1.8
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
    axs_kde = [
        [fig.add_subplot(gs[4, :2]), fig.add_subplot(gs[4, 3:5])],
        [fig.add_subplot(gs[5, :2]), fig.add_subplot(gs[5, 3:5])],
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

            sns.stripplot(
                data=df_agg[df_agg["marker"] == marker],
                x="stage",
                y=feature,
                hue="condition",
                order=["8hpf", "24hpf"],
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
                data=df_agg[df_agg["marker"] == marker],
                x="stage",
                y=feature,
                hue="condition",
                order=["8hpf", "24hpf"],
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

            ax[ax_idx].set_title(f"Mean", fontweight="bold")

    for feature, ax in zip(features_sum, axs_sum):
        for ax_idx, marker in enumerate(["gra", "dnd1"]):
            if feature in log_scale_features:
                ax[ax_idx].set_yscale("log")

            sns.stripplot(
                data=df_agg_sum[df_agg_sum["marker"] == marker],
                x="stage",
                y=feature,
                hue="condition",
                order=["8hpf", "24hpf"],
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
                data=df_agg_sum[df_agg_sum["marker"] == marker],
                x="stage",
                y=feature,
                hue="condition",
                order=["8hpf", "24hpf"],
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

            ax[ax_idx].set_title(f"Sum", fontweight="bold")

    for feature, ax in zip(features_std, axs_std):
        for ax_idx, marker in enumerate(["gra", "dnd1"]):
            if feature in log_scale_features:
                ax[ax_idx].set_yscale("log")

            sns.stripplot(
                data=df_agg_std[df_agg_std["marker"] == marker],
                x="stage",
                y=feature,
                hue="condition",
                order=["8hpf", "24hpf"],
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
                data=df_agg_std[df_agg_std["marker"] == marker],
                x="stage",
                y=feature,
                hue="condition",
                order=["8hpf", "24hpf"],
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
                err_kws=dict(linewidth=0.4, color="black"),
                markeredgewidth=1,
                palette="dark:black",
                zorder=5,
            )

            ax[ax_idx].set_title(f"Std", fontweight="bold")

    def _forward(x):
        return np.where(
            x < 0,
            (x + 1) * nuclear_fraction,
            nuclear_fraction + x * (1 - nuclear_fraction),
        )

    def _inverse(x):
        return np.where(
            x < nuclear_fraction,
            x / nuclear_fraction - 1,
            (x - nuclear_fraction) / (1 - nuclear_fraction),
        )

    for ax, dev_stage in zip(axs_kde, ["8hpf", "24hpf"]):
        for ax_idx, marker in enumerate(["gra", "dnd1"]):
            data_kde = df_agg[
                (df_agg["marker"] == marker) & (df_agg["stage"] == dev_stage)
            ]
            nuclear_fraction = float(data_kde["FranctionAlongRayToCellBoundary"].mean())

            print(
                f"Marker: {marker}, Stage: {dev_stage}, Nuclear Fraction: {nuclear_fraction:.2f}"
            )

            sns.kdeplot(
                data=data_kde,
                x="RelativeNuclearDistance",
                hue="condition",
                palette=color_palette,
                fill=True,
                ax=ax[ax_idx],
            )
            ax[ax_idx].set_xscale("function", functions=(_forward, _inverse))
            ax[ax_idx].set_xlim(-0.5, 1)
            tick_positions = [-0.5, 0, 0.5, 1]
            # tick_labels = ["Nucleus centroid", "Nucleus boundary", "Cell boundary"]
            ax[ax_idx].set_xlabel("Relative Nuclear Distance")
            ax[ax_idx].set_xticks(tick_positions)  # , labels=tick_labels)
            aspect_ratio = 0.05
            ax[ax_idx].set_aspect(aspect_ratio)
            ax[ax_idx].set_title(
                f"Relative Nuclear Distance of {marker.capitalize()} at {dev_stage}",
                fontweight="bold",
            )

    # Add a vertical dashed line between the two subplots
    mid_x = 0.5  # Horizontal midpoint of the figure (between the two subplots)
    line = mlines.Line2D(
        [mid_x, mid_x],
        [0.1, 0.9],  # x coords, y coords in figure space
        transform=fig.transFigure,
        color="gray",
        linestyle="--",
        linewidth=1,
    )
    fig.add_artist(line)
    plt.tight_layout()
    sns.despine()
    figure_dir = save_dir / "figures" / "per_cell_features"
    figure_dir.mkdir(exist_ok=True, parents=True)
    fig.savefig(
        figure_dir / f"plot_per_granule_profile_aggregate-{aggregation_function}.pdf",
        transparent=True,
    )
    plt.close(fig)


def plot_foldchange(df, save_dir: Path, color_palette="Set2"):
    """
    Plot the differences between conditions for each marker and stage
    """
    excluded_columns = {
        "uid",
        "condition",
        "stage",
        "marker",
        "GranuleIndex",
        "FeatureCategory",
    }
    feature_columns = [
        column
        for column in df.select_dtypes(include=[np.number]).columns
        if column not in excluded_columns and df[column].notna().any()
    ]

    if isinstance(color_palette, dict):
        palette = list(color_palette.values())
    else:
        palette = sns.color_palette(color_palette, n_colors=2)

    figure_dir = save_dir / "figures" / "fold_change"
    figure_dir.mkdir(exist_ok=True, parents=True)

    stage_order = [
        stage for stage in ["8hpf", "24hpf"] if stage in df["stage"].unique()
    ]
    if not stage_order:
        stage_order = list(df["stage"].dropna().unique())

    marker_order = [
        marker for marker in ["gra", "dnd1"] if marker in df["marker"].unique()
    ]
    if not marker_order:
        marker_order = list(df["marker"].dropna().unique())

    conditions = [
        condition
        for condition in ["ctrl", "kd"]
        if condition in df["condition"].unique()
    ]

    for stage in stage_order:
        for marker in marker_order:
            subset = df[(df["stage"] == stage) & (df["marker"] == marker)]
            subset = (
                subset.groupby(["uid", "condition"])[feature_columns]
                .mean()
                .reset_index()
            )
            if subset.empty:
                continue

            x_positions = np.arange(len(feature_columns))
            fold_changes = []
            p_values = []
            summary_lines = [f"{marker} at {stage}:"]
            second_col_start = 35
            summary_lines.append(
                f"   n_samples ctrl:{' ' * (second_col_start - 14)}{len(subset.loc[subset['condition'] == 'ctrl'])}"
            )
            summary_lines.append(
                f"   n_samples kd:{' ' * (second_col_start - 12)}{len(subset.loc[subset['condition'] == 'kd'])}"
            )

            def _clean_feature_values(series):
                values = series.dropna().to_numpy()
                if np.iscomplexobj(values):
                    values = np.real(values)
                values = np.asarray(values, dtype=float)
                return values[np.isfinite(values)]

            for feature in feature_columns:
                ctrl_values = _clean_feature_values(
                    subset.loc[subset["condition"] == "ctrl", feature]
                )
                kd_values = _clean_feature_values(
                    subset.loc[subset["condition"] == "kd", feature]
                )

                if len(ctrl_values) == 0 or len(kd_values) == 0:
                    fold_changes.append(np.nan)
                else:
                    ctrl_mean = float(np.mean(ctrl_values))
                    kd_mean = float(np.mean(kd_values))
                    fold_changes.append(
                        abs(float(kd_mean / ctrl_mean)) if ctrl_mean != 0 else np.nan
                    )

                if len(conditions) == 2:
                    if len(ctrl_values) < 2 or len(kd_values) < 2:
                        p_values.append(np.nan)
                    else:
                        _, _, p_value = statistical_analysis(
                            ctrl_values.tolist(), kd_values.tolist()
                        )
                        p_values.append(float(p_value))
                else:
                    p_values.append(np.nan)

                ratio_text = (
                    "nan" if np.isnan(fold_changes[-1]) else f"{fold_changes[-1]:.6g}"
                )
                p_value_text = (
                    "nan" if np.isnan(p_values[-1]) else f"{p_values[-1]:.6g}"
                )
                summary_lines.append(
                    f"   {feature}:{' ' * (second_col_start - len(str(feature)))}{ratio_text[:5]}\t {p_value_text}\t {get_stars(p_values[-1])}"
                )

            fig, (ax_mean, ax_pvalue) = plt.subplots(
                2,
                1,
                figsize=(max(12, len(feature_columns) * 0.5), 8),
                sharex=True,
            )

            ax_mean.bar(
                x_positions,
                fold_changes,
                color=palette[0],
                edgecolor="black",
                alpha=0.8,
            )
            ax_mean.set_ylabel("kd / ctrl")
            ax_mean.set_title(f"{stage} / {marker} - feature fold change")
            ax_mean.tick_params(axis="x", labelbottom=False)
            ax_mean.set_ylim(0, 10)
            ax_mean.axhline(1, color="black", linestyle="--", linewidth=1)

            ax_pvalue.bar(
                x_positions,
                p_values,
                color=palette[1] if len(palette) > 1 else palette[0],
                alpha=0.8,
                edgecolor="black",
            )
            ax_pvalue.axhline(0.05, color="black", linestyle="--", linewidth=1)
            ax_pvalue.set_ylabel("p-value")
            ax_pvalue.set_xlabel("Features")
            ax_pvalue.set_ylim(10e-10, 1)
            ax_pvalue.set_xticks(x_positions)
            ax_pvalue.set_xticklabels(feature_columns, rotation=45, ha="right")
            ax_pvalue.set_title(f"{stage} / {marker} - ctrl vs kd p-values")
            ax_pvalue.set_yscale("log")

            for x_pos, p_value in zip(x_positions, p_values):
                if np.isnan(p_value):
                    label_y = 1e-9
                else:
                    label_y = max(p_value * 1.4, 1e-9)
                ax_pvalue.text(
                    x_pos,
                    label_y,
                    get_stars(p_value),
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    fontweight="bold",
                )

            fig.tight_layout()
            sns.despine(fig=fig)
            fig.savefig(
                figure_dir / f"fold_change_{stage}_{marker}.pdf", bbox_inches="tight"
            )
            plt.close(fig)

            summary_path = figure_dir / f"fold_change_{stage}_{marker}.txt"
            summary_path.write_text("\n".join(summary_lines) + "\n")


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
    figure_dir = save_dir / "figures" / "individial_granule_data"
    figure_dir.mkdir(exist_ok=True, parents=True)
    fig.savefig(figure_dir / "ridgeplots_per_marker.pdf")
    plt.close(fig)


def cell_features(df, save_dir: Path, color_palette="Set2"):
    """
    Plot the cell features for each marker and stage
    """
    skip_columns = [
        "uid",
        "markers",
        "condition",
        "stage",
    ]
    df_cell = df.copy()

    conditions = df_cell["condition"].dropna().unique()
    stages = df_cell["stage"].dropna().unique()

    feature_columns = [
        column
        for column in df_cell.select_dtypes(include=[np.number]).columns
        if column not in skip_columns and df_cell[column].notna().any()
    ]

    for stage in stages:
        if not len(conditions) == 2:
            raise ValueError(
                f"Expected exactly 2 conditions for stage {stage}, but found {len(conditions)}: {conditions}"
            )
        feature_names = []
        feature_means = []
        feature_pvalues = []
        sample_numbers = []
        for feature in feature_columns:
            condition_values = []
            for condition in conditions:
                df_cond = df_cell[
                    (df_cell["stage"] == stage) & (df_cell["condition"] == condition)
                ]
                cond_values = df_cond[feature].dropna().to_numpy()
                condition_values.append(
                    {
                        "condition_name": condition,
                        "mean": np.mean(cond_values)
                        if len(cond_values) > 0
                        else np.nan,
                        "count": len(cond_values),
                        "values": cond_values.tolist(),
                    }
                )
            if not len(condition_values) == 2:
                raise ValueError(
                    f"Expected exactly 2 conditions for feature {feature} at stage {stage}, but found {len(condition_values)}: {condition_values}"
                )
            feature_names.append(feature)
            feature_means.append(
                abs(float(condition_values[1]["mean"] / condition_values[0]["mean"]))
            )
            feature_pvalues.append(
                get_stars(
                    statistical_analysis(
                        condition_values[0]["values"], condition_values[1]["values"]
                    )[2]
                )
            )
            sample_numbers.append(
                (condition_values[0]["count"], condition_values[1]["count"])
            )

        fig, ax = plt.subplots(
            2, 1, figsize=(max(12, len(feature_columns) * 0.5), 8), sharex=True
        )

        ax[0].bar(
            feature_names,
            feature_means,
            color="#B9B9B9",
            alpha=0.8,
            edgecolor="black",
        )

        ax[0].set_ylabel(f"{conditions[1]} / {conditions[0]}")
        ax[0].set_title(
            f"{stage} - cell feature fold change", fontsize=11, fontweight="bold"
        )
        ax[0].set_ylim(0, 10)
        ax[0].tick_params(axis="x", labelbottom=False)
        ax[1].bar(
            feature_names,
            feature_pvalues,
            color="#B9B9B9",
            alpha=0.8,
            edgecolor="black",
        )
        ax[1].axhline(0.05, color="black", linestyle="--", linewidth=1)
        ax[1].set_ylabel("p-value")
        ax[1].set_xlabel("Features")
        ax[1].set_ylim(10e-10, 1)
        ax[1].set_xticks(range(len(feature_names)))
        ax[1].set_xticklabels(feature_names, rotation=45, ha="right")
        ax[1].set_title(
            f"{stage} - ctrl vs kd p-values", fontsize=11, fontweight="bold"
        )
        ax[1].set_yscale("log")

        fig.tight_layout()
        save_path = (
            save_dir
            / "figures"
            / "cell_features"
            / f"cell_feature_fold_change_{stage}.pdf"
        )
        save_path.parent.mkdir(exist_ok=True, parents=True)
        fig.savefig(save_path)
        plt.close(fig)

        with open(save_path.parent / f"cell_feature_fold_change_{stage}.txt", "w") as f:
            f.write(f"{stage} cell feature fold change:\n")
            for name, mean, pval, (n1, n2) in zip(
                feature_names, feature_means, feature_pvalues, sample_numbers
            ):
                f.write(f"   {name}:   {mean:.6g}\t {pval}\t (n={n1}, n={n2})\n")
