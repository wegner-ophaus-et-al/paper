import sys

sys.path.insert(0, "/Users/julian/Documents/General Science/Programming/py/packages")

from gamgee.segmenter import SegmentationModel

from pathlib import Path
from cell import Cell
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from tqdm import tqdm
import pandas as pd
from utils import calculate_area_ratio, statistics

data_root = Path(__file__).parent / "data"

color_palette = {"DMSO": "#888888", "CHX": "#66B2DD", "PatA": "#E6A925"}

brightness_threshold = 11000
# brightness_threshold = np.percentile(df_int["value"], 50)

txt_results = []
cells = []
for lsm_path in data_root.rglob("*.lsm"):
    if not lsm_path.name.startswith("."):
        cells.append(Cell(lsm_path.parent.parent))

segment_new = False

sm_granule = None
if segment_new:
    sm_granule = SegmentationModel(
        "/Users/icb_remote/Documents/JW/py/packages/gamgee/models/msam/sam_granules_refined_up3_34917658",
        model_type="vit_l_lm",
    )


ncols = 3
nrows = len(cells)
width_size = 1.63
height_size = 2
fig, axes = plt.subplots(
    nrows=len(cells), ncols=3, figsize=(width_size * ncols, height_size * nrows)
)
fig_segmentation, axes_segmentation = plt.subplots(
    nrows=len(cells), ncols=6, figsize=(width_size * 6, height_size * nrows)
)


data = []
with tqdm(total=len(cells), desc="Processing cells") as pbar:
    for ax, ax_seg, cell in zip(axes, axes_segmentation, cells):
        cell.read_images()  # Load images from disk
        if segment_new:
            cell.segment(sm_granule)  # Uncomment this line to perform segmentation
            cell.write_images()  # Uncomment this line if segmentations are performed
        cell.clean_segmentations()
        cell.plot_images(ax)
        cell.generate_various_masks()
        cell.plot_segmentations(ax_seg)
        data.extend(cell.eval())
        pbar.update(1)

plt.tight_layout()
fig.savefig(data_root / "contactsheet_segmentation_results.pdf", dpi=300)
plt.close(fig)
plt.tight_layout()
fig_segmentation.savefig(data_root / "contactsheet_generated_masks.pdf", dpi=100)
plt.close(fig_segmentation)

df_prime = pd.DataFrame(data)
df_prime.to_csv(data_root / "measurements_all_cell.csv", index=False)
if brightness_threshold is None:
    df = df_prime.copy()
else:
    df_int = df_prime[
        (df_prime["segmentation_name"] == "granule")
        & (df_prime["image_name"] == "granule")
        & (df_prime["measurement_name"] == "mean")
    ].copy()
    bright_cells = df_int[df_int["value"] > brightness_threshold]["uid"].unique()
    df = df_prime[df_prime["uid"].isin(bright_cells)].copy()
    df.to_csv(data_root / "measurements_bright_cells.csv", index=False)

unique_measurements = df["measurement_name"].unique()
number_unique_measurements = len(unique_measurements)


def generate_figures(df):
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
    unique_measurements_fig = df["measurement_name"].unique()
    number_unique_measurements_fig = len(unique_measurements_fig)
    ncols = 8
    nrows = number_unique_measurements_fig

    grixel_width = 1.63
    grixel_height = 2

    fig_data_fig = plt.figure(figsize=(grixel_width * ncols, grixel_height * nrows))
    gs = GridSpec(nrows, ncols)
    axs_fig = []
    for i in range(number_unique_measurements_fig):
        axs_fig.append(
            [
                [
                    fig_data_fig.add_subplot(gs[i, :3]),
                    fig_data_fig.add_subplot(gs[i, 3]),
                ],  # actin
                [
                    fig_data_fig.add_subplot(gs[i, 4:7]),
                    fig_data_fig.add_subplot(gs[i, 7]),
                ],  # vasa
            ]
        )
    return fig_data_fig, axs_fig


for channel in df["image_name"].unique():
    sns.set_style("ticks")
    fig_data, axs = generate_figures(df)
    df_channel = df[df["image_name"] == channel]
    for unique_measurement, ax in zip(unique_measurements, axs):
        df_measurement = df_channel[
            df_channel["measurement_name"] == unique_measurement
        ]

        txt_results.append(f"Results for {channel}-{unique_measurement}:\n")

        # actin plots
        sns.stripplot(
            data=df_measurement[df_measurement["rna_type"] == "actin"],
            y="value",
            x="segmentation_name",
            hue="condition",
            palette=color_palette,
            dodge=0.5,
            alpha=0.4,
            size=3,
            ax=ax[0][0],
        )
        sns.pointplot(
            data=df_measurement[df_measurement["rna_type"] == "actin"],
            y="value",
            x="segmentation_name",
            hue="condition",
            dodge=0.6,
            errorbar="sd",
            estimator="median",
            capsize=0.075,
            linestyle="none",
            markersize=10,
            marker="_",
            err_kws=dict(linewidth=0.4, color="black"),
            markeredgewidth=1,
            palette="dark:black",
            zorder=5,
            ax=ax[0][0],
        )
        ax[0][0].set_title(f"{channel} - {unique_measurement} - actin")
        ax[0][0].set_ylabel(f"{channel} signal in masked area (a.u.)")
        ax[0][0].set_xlabel("Masked area")

        ratio_col_a = "nucleo_cytoplasm"
        ratio_col_b = "granule_periphery"
        df_ratio_actin = calculate_area_ratio(
            df[df["rna_type"] == "actin"],
            channel,
            ratio_col_a,
            ratio_col_b,
            measurement=unique_measurement,
        )
        txt_results.extend(
            statistics(df_ratio_actin, "actin", f"{ratio_col_a} / {ratio_col_b}")
        )

        sns.stripplot(
            df_ratio_actin.reset_index(),
            y="ratio",
            x="condition",
            hue="condition",
            palette=color_palette,
            alpha=0.4,
            size=3,
            ax=ax[0][1],
        )
        sns.pointplot(
            data=df_ratio_actin.reset_index(),
            x="condition",
            y="ratio",
            hue="condition",
            # dodge=0.4,
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
            ax=ax[0][1],
        )

        ax[0][1].set_title(f"{ratio_col_a} / {ratio_col_b} ratio")
        ax[0][1].set_ylabel(f"Ratio of {channel} signal {ratio_col_a} / {ratio_col_b}")
        ax[0][1].set_xlabel("Condition")

        # vasa plots

        sns.stripplot(
            data=df_measurement[df_measurement["rna_type"] == "vasa"],
            y="value",
            x="segmentation_name",
            hue="condition",
            dodge=True,
            palette=color_palette,
            alpha=0.4,
            size=3,
            ax=ax[1][0],
        )
        sns.pointplot(
            data=df_measurement[df_measurement["rna_type"] == "vasa"],
            y="value",
            x="segmentation_name",
            hue="condition",
            dodge=0.6,
            errorbar="sd",
            estimator="median",
            capsize=0.075,
            linestyle="none",
            markersize=10,
            marker="_",
            err_kws=dict(linewidth=0.4, color="black"),
            markeredgewidth=1,
            palette="dark:black",
            zorder=5,
            ax=ax[1][0],
        )
        ax[1][0].set_title(f"{channel} - {unique_measurement} - vasa")
        ax[1][0].set_ylabel(f"{channel} signal in masked area (a.u.)")
        ax[1][0].set_xlabel("Masked area")

        df_ratio_vasa = calculate_area_ratio(
            df[df["rna_type"] == "vasa"],
            channel,
            ratio_col_a,
            ratio_col_b,
            measurement=unique_measurement,
        )

        txt_results.extend(
            statistics(df_ratio_vasa, "vasa", f"{ratio_col_a} / {ratio_col_b}")
        )

        sns.stripplot(
            df_ratio_vasa.reset_index(),
            y="ratio",
            x="condition",
            hue="condition",
            palette=color_palette,
            alpha=0.4,
            size=3,
            ax=ax[1][1],
        )
        sns.pointplot(
            data=df_ratio_vasa.reset_index(),
            x="condition",
            y="ratio",
            hue="condition",
            # dodge=0.4,
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
            ax=ax[1][1],
        )
        ax[1][1].set_title(f"{ratio_col_a} / {ratio_col_b} ratio")
        ax[1][1].set_ylabel(f"Ratio of {channel} signal {ratio_col_a} / {ratio_col_b}")
        ax[1][1].set_xlabel("Condition")
        # ax[1][1].set_ylim(0, 2)

        for ax_row in ax:
            for ax_single in ax_row:
                # remove legend
                # ax_single.legend_.remove()
                pass

        for rna_tpe in ["actin", "vasa"]:
            for ratio_comb in [
                ("granule", "granule_periphery"),
                ("cytoplasm", "granule_periphery"),
                ("nucleo_cytoplasm", "granule"),
            ]:
                txt_results.extend(
                    statistics(
                        calculate_area_ratio(
                            df[df["rna_type"] == rna_tpe],
                            channel,
                            ratio_comb[0],
                            ratio_comb[1],
                            measurement=unique_measurement,
                        ),
                        rna_tpe,
                        f"{ratio_comb[0]} / {ratio_comb[1]}",
                    )
                )

    sns.despine()
    plt.tight_layout()
    fig_data_path = data_root / "figures"
    fig_data_path.mkdir(exist_ok=True)
    fig_data.savefig(fig_data_path / f"{channel}_measurements.pdf", transparent=True)

with open(data_root / "statistical_results.txt", "w") as f:
    f.write("\n".join(txt_results))
