import tifffile as tiff
from thecell import TheCell, mh
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
from utils import print_nested, statistical_analysis, get_stars
from multiprocessing import Pool
import pandas as pd
from tqdm import tqdm
import seaborn as sns


def compute_sample(sample_dir):

    original_raw_img = tiff.imread(next((sample_dir / "original_raw").glob("*.tif")))

    img = np.sum(original_raw_img, axis=0)

    sum_images = {}
    if "2025-03-13" in sample_dir.name:
        sum_images["dnd"] = img[0]
        sum_images["gra"] = img[2]
        sum_images["cell"] = img[1]
        sum_images["nucleus"] = img[1]

    else:
        sum_images["cell"] = img[0]
        sum_images["nucleus"] = img[0]
        sum_images["dnd"] = img[1]
        sum_images["gra"] = img[2]

    cell_obj = TheCell(
        sample_dir,
        conditions=["full_mix", "no_tdrd7a"],
        model_handler=mh,
        granuleA="dnd",
        granuleB="gra",
    )

    cell_obj.read_segmentations()
    # Replace the raw_image in each marker with the corresponding image from the sum projection
    for marker_name, marker in cell_obj.markers.items():
        cell_obj.markers[marker_name].raw_image = sum_images[marker_name]

    cell_obj.compute_features()
    fig, ax = plt.subplots(1, 3, figsize=(12, 4))
    ax[0].imshow(sum_images["gra"], cmap="gray")
    ax[0].set_title("GRA")
    ax[1].imshow(sum_images["dnd"], cmap="gray")
    ax[1].set_title("DND")
    ax[2].imshow(sum_images["cell"], cmap="gray")
    ax[2].set_title("Cell/Nucleus")

    ax[0].contour(cell_obj.markers["gra"].segmentation, colors="cyan", linewidths=0.5)
    ax[1].contour(cell_obj.markers["dnd"].segmentation, colors="yellow", linewidths=0.5)
    ax[2].contour(
        cell_obj.markers["cell"].segmentation, colors="magenta", linewidths=0.5
    )

    fig.savefig(
        Path(__file__).parent / "proto_out" / f"{cell_obj.uid}_sum_projection.png",
        dpi=300,
    )
    return cell_obj


def main():
    root_path = Path(
        "/Users/icb_remote/Documents/JW/py/data/iPGC_24hpf_tdrd7a-modulation/fm_nt_data"
    )
    sample_dirs = [
        d for d in root_path.iterdir() if d.is_dir() and "figures" not in d.name
    ]
    cell_objects = []
    pbar = tqdm(total=len(sample_dirs), desc="Processing Samples", unit="sample")
    with Pool(processes=5) as pool:
        cell_objects = list(pool.map(compute_sample, sample_dirs))
        pbar.update(1)

    cell_feature_dicts = list(map(lambda cell: cell.get_cell_features(), cell_objects))

    df = pd.DataFrame(cell_feature_dicts)
    df.to_csv(Path(__file__).parent / "cell_features.csv", index=False)

    stat_outpuit_lines = []

    color_palette = {
        "full_mix": "#878787",
        "no_tdrd7a": "#8a38a6",
    }
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

    ncols = 8
    nrows = 6
    subfig_cw = 1.4
    subfig_ch = 1.6

    fig = plt.figure(figsize=(ncols * subfig_cw, subfig_ch * nrows))

    gs = gridspec.GridSpec(ncols=ncols, nrows=nrows)

    axs = []

    for i in range(nrows):
        axs.append(
            [
                [fig.add_subplot(gs[i, 0]), fig.add_subplot(gs[i, 4])],
                [fig.add_subplot(gs[i, 1]), fig.add_subplot(gs[i, 5])],
                [fig.add_subplot(gs[i, 2]), fig.add_subplot(gs[i, 6])],
                [fig.add_subplot(gs[i, 3]), fig.add_subplot(gs[i, 7])],
            ]
        )

    for ax_row, stat in zip(axs, ["Sum", "Mean", "Std", "Skew", "Kurtosis", "Cv"]):
        for ax_region, region in zip(
            ax_row, ["Granules", "Nucleus", "Cytoplasm", "Cell"]
        ):
            for ax, marker_name in zip(ax_region, ["gra", "dnd"]):
                feature_name = f"{marker_name}_{stat}{region}"

                sns.stripplot(
                    data=df,
                    # x="stage",
                    y=feature_name,
                    hue="condition",
                    # order=["8hpf", "24hpf"],
                    hue_order=list(color_palette.keys())
                    if isinstance(color_palette, dict)
                    else None,
                    palette=color_palette,
                    dodge=0.4,
                    alpha=0.4,
                    size=2,
                    ax=ax,
                    jitter=0.3,
                )

                sns.pointplot(
                    data=df,
                    # x="stage",
                    y=feature_name,
                    hue="condition",
                    # order=["8hpf", "24hpf"],
                    hue_order=list(color_palette.keys())
                    if isinstance(color_palette, dict)
                    else None,
                    dodge=0.4,
                    ax=ax,
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

                ax.set_title(f"{marker_name.upper()} {stat} {region}")

                fm_values = (
                    df[df["condition"] == "full_mix"][feature_name].dropna().to_list()
                )
                nt_values = (
                    df[df["condition"] == "no_tdrd7a"][feature_name].dropna().to_list()
                )

                test_name, _, p_value = statistical_analysis(fm_values, nt_values)
                star_significance = get_stars(p_value)
                printed_p_val = (
                    # str(round(p_value, 4)) if p_value >= 0.0001 else "<0.0001"
                    p_value
                )
                stat_outpuit_lines.append(
                    f"{feature_name}:{(30 - len(feature_name)) * ' '}{printed_p_val}\t{star_significance}\t{test_name}\n"
                )

    fig.tight_layout()
    sns.despine()
    fig.savefig(
        root_path / "figures" / "cell_features" / "sum_projection_cell_features.pdf"
    )
    with open(
        root_path / "figures" / "cell_features" / "statistical_analysis.txt", "w"
    ) as f:
        for line in stat_outpuit_lines:
            f.write(line)


if __name__ == "__main__":
    main()
