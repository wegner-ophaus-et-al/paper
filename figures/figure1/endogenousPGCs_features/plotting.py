import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
from pathlib import Path


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


def plot_individial_granule_profile(df, save_dir: Path, color_paltte="Set2"):
    """
    Plot the granule profile for all cells in the dataframe
    """

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
            palette=color_paltte,
            ax=ax,
            split=True,
            inner="quart",
        )
    fig.tight_layout()
    sns.despine()
    figure_dir = save_dir / "figures"
    figure_dir.mkdir(exist_ok=True)
    fig.savefig(figure_dir / "Individual.pdf")
