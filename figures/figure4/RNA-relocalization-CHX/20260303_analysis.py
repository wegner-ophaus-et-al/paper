from skimage import measure
from scipy.ndimage import distance_transform_edt
from scipy import stats
import matplotlib.pyplot as plt
from merge_utils import merge_images_with_cmaps
from colormaps import cyan, magenta, yellow, green, red
import tifffile as tiff
import numpy as np
from pathlib import Path
import seaborn as sns
import pandas as pd
from multiprocessing import Pool


def parametric(data_set1, data_set2):
    set1_normal = stats.normaltest(data_set1).pvalue > 0.05
    set2_normal = stats.normaltest(data_set2).pvalue > 0.05

    return min(set1_normal, set2_normal)


def statistical_analysis(data_set1, data_set2):
    if parametric(data_set1, data_set2):
        t_statistic, p_value = stats.ttest_ind(data_set1, data_set2)
        test_type = "t-test"
    else:
        t_statistic, p_value = stats.mannwhitneyu(data_set1, data_set2)
        test_type = "Mann-Whitney U test"
    return test_type, t_statistic, p_value


def generate_granule_periphery_mask(seg_granules, distance_in=3, distance_out=5):
    """
    Generate a mask for the periphery region around granules based on distance transform.
    :param seg_granules: Segmented granules image
    :param distance_in: Distance in pixels from the edge towards the center (inner boundary)
    :param distance_out: Distance in pixels from the edge towards the outside (outer boundary)
    :return: Binary mask of the periphery region
    """
    # Create a binary mask of granules
    granule_mask = seg_granules > 0

    # Compute distance transform INSIDE granules
    distance_inside = distance_transform_edt(granule_mask)

    # Compute distance transform OUTSIDE granules
    distance_outside = distance_transform_edt(~granule_mask)

    # Create periphery mask: inside boundary + outside boundary
    periphery_mask = (distance_inside <= distance_in) & (
        distance_outside <= distance_out
    )

    return periphery_mask


def process_sample(
    sample_dir,
    ax=None,
    percentile_vasa=99.0,
    percentile_rna=99.0,
    save_images=False,
    **kwargs,
):
    """
    Process a single sample directory to extract intensity and morphological features.
    :param sample_dir:
    :param ax:
    :param percentile_vasa:
    :param percentile_rna:
    :param save_images:
    :param kwargs: distance_in=int, distance_out=int for periphery mask number of pixels of periphery to the inside and outside
    :return:
    """
    print(f"Processing sample: {sample_dir.name}")

    uid, condition, name = sample_dir.name.split("__")
    name_splits = name.split("_")

    date = name_splits[0]  # Fist fraction is the date YYYY-MM-DD
    cell_name = name_splits[-1]  # Last fraction is the cell name e.g. e2c1
    repeat_name = name_splits[
        1:-1
    ]  # Middle fractions the name of the repeat e.g. DMSA_b2

    condition_translate = {"ctrl": "DMSO Control", "chx": "CHX Treated"}

    img_granules = tiff.imread(sample_dir / "imgs" / "vasa.tif")
    img_nanos_rna = tiff.imread(sample_dir / "imgs" / "nanos.tif")
    img_tdrd7_rna = tiff.imread(sample_dir / "imgs" / "tdrd7a.tif")

    seg_cell = tiff.imread(sample_dir / "segmentations" / "cell.tif")
    seg_granules = tiff.imread(sample_dir / "segmentations" / "vasa.tif")
    seg_nanos_rna = tiff.imread(sample_dir / "segmentations" / "nanos.tif")
    seg_tdrd7_rna = tiff.imread(sample_dir / "segmentations" / "tdrd7a.tif")

    seg_granules = measure.label(
        seg_granules
    )  # .astype(bool) # Ensure granule segmentation is binary
    # Generate periphery mask around granules
    periphery_mask = generate_granule_periphery_mask(seg_granules, **kwargs)
    # Mask for the granule centers (granule area excluding periphery)
    granule_center_mask = (seg_granules > 0) & (~periphery_mask)

    # Mask for cytoplasm (cell excluding granules and periphery)
    cytoplasm_mask = (seg_cell > 0) & (seg_granules == 0) & (~periphery_mask)

    results_dict = dict()
    results_list = []

    # Initialze resutls dict with metadata
    results_dict["UID"] = uid
    results_dict["Condition"] = condition_translate.get(condition, condition)
    results_dict["Date"] = date
    results_dict["CellName"] = cell_name
    results_dict["RepeatName"] = "_".join(repeat_name)

    # Get granule count and sizes
    granule_props = measure.regionprops(seg_granules)
    granule_sizes = [prop.area for prop in granule_props]
    granule_mean_size = np.mean(granule_sizes) if granule_sizes else 0
    results_dict["GranuleCount"] = len(granule_sizes)
    results_dict["GranuleSizes"] = granule_sizes
    results_dict["MeanGranuleSize"] = granule_mean_size

    # Calculate Vasa metrics
    results_dict["MeanVasaInCell"] = np.mean(img_granules[seg_cell > 0])
    results_dict["TotalVasaInCell"] = np.sum(img_granules[seg_cell > 0])

    # Calculate mean and total intensities in different regions for each RNA
    for rna_name, rna_img, rna_seg in zip(
        ["nanos", "tdrd7"],
        [img_nanos_rna, img_tdrd7_rna],
        [seg_nanos_rna, seg_tdrd7_rna],
    ):
        temp_dict = results_dict.copy()
        temp_dict["RNAType"] = rna_name
        temp_dict["MeanRNAInGranules"] = np.mean(rna_img[seg_granules > 0])
        temp_dict["MeanRNAInPeriphery"] = np.mean(rna_img[periphery_mask])
        temp_dict["MeanRNAInCytoplasm"] = np.mean(rna_img[cytoplasm_mask])
        temp_dict["MeanRNAInGranuleCenters"] = np.mean(rna_img[granule_center_mask])

        temp_dict["TotalRNAInGranules"] = np.sum(rna_img[seg_granules > 0])
        temp_dict["TotalRNAInPeriphery"] = np.sum(rna_img[periphery_mask])
        temp_dict["TotalRNAInCytoplasm"] = np.sum(rna_img[cytoplasm_mask])
        temp_dict["TotalRNAInGranuleCenters"] = np.sum(rna_img[granule_center_mask])

        # Get RNA puncta counts
        rna_props = measure.regionprops(rna_seg)
        temp_dict["RNAPunctaCount"] = len(rna_props)

        # Get area of overlap between RNA puncta with granules
        rna_granule_overlap = (seg_granules > 0) & (rna_seg > 0)
        temp_dict["RNAGranuleOverlapArea"] = np.sum(rna_granule_overlap)
        temp_dict["RNAGranuleTotalArea"] = np.sum(rna_seg > 0)

        results_list.append(temp_dict)

    # Visualize images, segmentations, and masks
    if ax is not None:
        if len(ax) != 6:
            raise ValueError("ax must have 6 subplots.")
        # Write sample name vertically on the left side
        ax[0].text(
            -0.1,
            0.5,
            f"{uid}_{condition}",
            va="center",
            ha="right",
            rotation=90,
            transform=ax[0].transAxes,
            fontsize=6,
        )
        # Granules image
        ax[0].imshow(
            img_granules,
            cmap=cyan,
            vmin=0,
            vmax=np.percentile(img_granules, percentile_vasa),
        )  # Granules image

        # Nanos RNA image with granule overlay
        ax[1].imshow(
            img_nanos_rna,
            cmap=yellow,
            vmin=0,
            vmax=np.percentile(img_nanos_rna, percentile_rna),
        )  # Nanos RNA image
        ax[1].contour(
            seg_nanos_rna > 0, colors="red", linewidths=0.5
        )  # Overlay granule segmentation
        # Tdrd7 RNA image with granule overlay
        ax[2].imshow(
            img_tdrd7_rna,
            cmap=magenta,
            vmin=0,
            vmax=np.percentile(img_tdrd7_rna, percentile_rna),
        )  # Tdrd7 RNA image
        ax[2].contour(
            seg_tdrd7_rna > 0, colors="red", linewidths=0.5
        )  # Overlay granule segmentation
        # Segmented granules with periphery and cytoplasm masks
        ax[3].imshow(periphery_mask > 0, cmap="gray")  # Granule segmentation
        # ax[3].contour(periphery_mask, colors='blue', linewidths=0.5) # Overlay periphery mask
        ax[3].contour(
            cytoplasm_mask, colors="green", linewidths=0.5
        )  # Overlay cytoplasm mask
        # Nanos RNA segmentation with granule overlay
        ax[4].imshow(seg_nanos_rna, cmap="gray")  # Nanos RNA segmentation
        ax[4].contour(
            seg_granules > 0, colors="red", linewidths=0.5
        )  # Overlay granule segmentation
        # Tdrd7 RNA segmentation with granule overlay
        ax[5].imshow(seg_tdrd7_rna, cmap="gray")  # Tdrd7 RNA segmentation
        ax[5].contour(
            seg_granules > 0, colors="red", linewidths=0.5
        )  # Overlay granule segmentation

        for a in ax:
            a.axis("off")

    if save_images:
        fig_save, ax_save = plt.subplots(1, 4, figsize=(9, 3))

        # Granules image with dashed cell contour and semi-translucent periphery mask
        ax_save[0].imshow(
            img_granules,
            cmap=cyan,
            vmin=0,
            vmax=np.percentile(img_granules, percentile_vasa),
        )  # Granules image
        ax_save[0].imshow(
            periphery_mask > 0, cmap=yellow, alpha=0.3
        )  # Overlay periphery mask with transparency
        ax_save[0].contour(
            seg_cell > 0, colors="white", linewidths=0.5, linestyles="dashed"
        )  # Overlay cell contour
        ax_save[0].set_title("Granules Image", fontsize=12)

        # Nanos RNA image with segmentation contour
        ax_save[1].imshow(
            img_nanos_rna,
            cmap=yellow,
            vmin=0,
            vmax=np.percentile(img_nanos_rna, percentile_rna),
        )  # Nanos RNA image
        ax_save[1].contour(
            seg_nanos_rna > 0, colors="white", linewidths=0.5, linestyles="dashed"
        )  # Overlay granule segmentation
        ax_save[1].set_title("Nanos RNA Image", fontsize=12)

        # Tdrd7 RNA image with segmentation contour
        ax_save[2].imshow(
            img_tdrd7_rna,
            cmap=magenta,
            vmin=0,
            vmax=np.percentile(img_tdrd7_rna, percentile_rna),
        )  # Tdrd7 RNA image
        ax_save[2].contour(
            seg_tdrd7_rna > 0, colors="white", linewidths=0.5, linestyles="dashed"
        )  # Overlay granule segmentation
        ax_save[2].set_title("Tdrd7 RNA Image", fontsize=12)

        # Merge of the two RNA images with the granule segmentation contour
        merged_rna = merge_images_with_cmaps(
            img_tdrd7_rna,
            img_nanos_rna,
            magenta,
            yellow,
            vmin1=0,
            vmax1=np.percentile(img_tdrd7_rna, percentile_rna),
            vmin2=0,
            vmax2=np.percentile(img_nanos_rna, percentile_rna),
        )
        ax_save[3].imshow(merged_rna)
        # ax_save[3].imshow(np.zeros_like(img_nanos_rna), cmap='gray', vmin=0, vmax=1) # Blank background
        # ax_save[3].imshow(img_nanos_rna, cmap=yellow, vmin=0, vmax=np.percentile(img_nanos_rna, percentile_rna),
        #                   alpha=img_nanos_rna/img_nanos_rna.max())
        # ax_save[3].imshow(img_tdrd7_rna, cmap=magenta, vmin=0, vmax=np.percentile(img_tdrd7_rna, percentile_rna),
        #                   alpha=img_tdrd7_rna/img_tdrd7_rna.max())
        ax_save[3].contour(
            seg_granules > 0, colors="#00FFFF", linewidths=0.5, linestyles="dashed"
        )  # Overlay granule segmentation
        ax_save[3].set_title("Merged RNA Images", fontsize=12)

        for a in ax_save:
            a.axis("off")
        plt.tight_layout()
        out_path = sample_dir / "output"
        out_path.mkdir(exist_ok=True, parents=True)
        fig_save.savefig(out_path / "summary_images.pdf", dpi=600)
        plt.close(fig_save)

    return results_list


root = Path("/Users/julian/local_files/chx_kim/data")
list_sample_dirs = sorted(
    [d for d in root.iterdir() if d.is_dir() and not d.name.startswith(".")]
)

fig, axes = plt.subplots(
    len(list_sample_dirs), 6, figsize=(18, 3 * len(list_sample_dirs))
)
result_collector = []
for i, sample_directory in enumerate(list_sample_dirs):
    axe = axes[i] if len(list_sample_dirs) > 1 else axes
    # Add header labels
    if i == 0:
        axe[0].set_title("Granules Image", fontsize=12)
        axe[1].set_title("Nanos RNA Image", fontsize=12)
        axe[2].set_title("Tdrd7 RNA Image", fontsize=12)
        axe[3].set_title("Granule Segmentation", fontsize=12)
        axe[4].set_title("Nanos RNA Segmentation", fontsize=12)
        axe[5].set_title("Tdrd7 RNA Segmentation", fontsize=12)
    result_collector.extend(
        process_sample(
            sample_directory,
            ax=axe,
            percentile_vasa=99.9,
            percentile_rna=99.4,
            distance_in=3,
            distance_out=5,
            save_images=False,
        )
    )

plt.tight_layout()
fig.savefig(root / "summary_analysis.pdf")

df = pd.DataFrame(result_collector)

df["Condition_RNA"] = df["Condition"] + "_" + df["RNAType"]

df_per_sample = (
    df.groupby(["UID", "Condition", "Date", "CellName", "RepeatName", "RNAType"])
    .agg(
        {
            "MeanRNAInGranules": "mean",
            "MeanRNAInPeriphery": "mean",
            "MeanRNAInCytoplasm": "mean",
            "MeanRNAInGranuleCenters": "mean",
            "TotalRNAInGranules": "mean",
            "TotalRNAInPeriphery": "mean",
            "TotalRNAInCytoplasm": "mean",
            "TotalRNAInGranuleCenters": "mean",
            "RNAPunctaCount": "mean",
            "RNAGranuleOverlapArea": "mean",
            "MeanGranuleSize": "mean",
            "GranuleCount": "mean",
        }
    )
    .reset_index()
)

# Calculate Various Ratios for mean intensities
df["MeanRNAInGranules_Cytoplasm_Ratio"] = (
    df["MeanRNAInGranules"] / df["MeanRNAInCytoplasm"]
)
df["MeanRNA_Cytoplasm_Granules_Ratio"] = (
    df["MeanRNAInCytoplasm"] / df["MeanRNAInGranules"]
)
df["MeanRNAInGranuleCenters_Cytoplasm_Ratio"] = (
    df["MeanRNAInGranuleCenters"] / df["MeanRNAInCytoplasm"]
)
df["MeanRNAInPeriphery_Cytoplasm_Ratio"] = (
    df["MeanRNAInPeriphery"] / df["MeanRNAInCytoplasm"]
)
df["MeanRNAInPeriphery_GranuleCenters_Ratio"] = (
    df["MeanRNAInPeriphery"] / df["MeanRNAInGranuleCenters"]
)
# Calculate Various Ratios for total intensities
df["TotalRNAInGranules_Cytoplasm_Ratio"] = (
    df["TotalRNAInGranules"] / df["TotalRNAInCytoplasm"]
)
df["TotalRNAInGranuleCenters_Cytoplasm_Ratio"] = (
    df["TotalRNAInGranuleCenters"] / df["TotalRNAInCytoplasm"]
)
df["TotalRNAInPeriphery_Cytoplasm_Ratio"] = (
    df["TotalRNAInPeriphery"] / df["TotalRNAInCytoplasm"]
)
df["TotalRNAInPeriphery_GranuleCenters_Ratio"] = (
    df["TotalRNAInPeriphery"] / df["TotalRNAInGranuleCenters"]
)

# Calculate Ratio of Segmentation overlap
df["RNAGranuleOverlap_Fraction"] = (
    df["RNAGranuleOverlapArea"] / df["RNAGranuleTotalArea"]
)

color_palette = {
    "DMSO Control_nanos": "#C8AF55",
    "DMSO Control_tdrd7": "#B94B4B",
    "CHX Treated_nanos": "#EDC001",
    "CHX Treated_tdrd7": "#F06F6F",
}


def get_significance_symbol(p_value):
    if p_value < 0.001:
        return "***"
    elif p_value < 0.01:
        return "**"
    elif p_value < 0.05:
        return "*"
    else:
        return "ns"


def get_signifcances(the_df, feature, rna_type):
    """
    Sort the dataset for the marker and perform a statistical test, comparing the control to the chx trewatment
    """
    sub_df = the_df.copy()
    sub_df = sub_df.query("RNAType == @rna_type")

    ctrl_vals = sub_df.query("Condition == 'DMSO Control'")[feature].values
    exp_vals = sub_df.query("Condition == 'CHX Treated'")[feature].values
    test_type, t_statistic, p_value = statistical_analysis(ctrl_vals, exp_vals)

    return {
        "feature": feature,
        "rna_type": rna_type,
        "test_type": test_type,
        "t_statistic": t_statistic,
        "p_value": p_value,
        "significance_symbol": get_significance_symbol(p_value),
    }


significance_testing = []

sns.set_style("ticks")
fig2, ax2 = plt.subplots(2, 4, figsize=(15, 10))
sns.violinplot(
    data=df,
    x="RNAType",
    y="MeanRNAInGranules_Cytoplasm_Ratio",
    hue="Condition_RNA",
    ax=ax2[0, 0],
    inner="quart",
    palette=color_palette,
)
sns.swarmplot(
    data=df,
    x="RNAType",
    y="MeanRNAInGranules_Cytoplasm_Ratio",
    hue="Condition_RNA",
    ax=ax2[0, 0],
    color="k",
    alpha=0.5,
    dodge=True,
)
ax2[0, 0].set_title("Mean Intensity (whole) Granules/Cytoplasm Ratio")

sns.violinplot(
    data=df,
    x="RNAType",
    y="MeanRNAInGranuleCenters_Cytoplasm_Ratio",
    hue="Condition_RNA",
    ax=ax2[0, 1],
    inner="quart",
    palette=color_palette,
)
sns.swarmplot(
    data=df,
    x="RNAType",
    y="MeanRNAInGranuleCenters_Cytoplasm_Ratio",
    hue="Condition_RNA",
    ax=ax2[0, 1],
    color="k",
    alpha=0.5,
    dodge=True,
)
ax2[0, 1].set_title("Mean Intensity (centers) Granules/Cytoplasm Ratio")


sns.violinplot(
    data=df,
    x="RNAType",
    y="MeanRNAInPeriphery_Cytoplasm_Ratio",
    hue="Condition_RNA",
    ax=ax2[0, 2],
    inner="quart",
    palette=color_palette,
)
sns.swarmplot(
    data=df,
    x="RNAType",
    y="MeanRNAInPeriphery_Cytoplasm_Ratio",
    hue="Condition_RNA",
    ax=ax2[0, 2],
    color="k",
    alpha=0.5,
    dodge=True,
)
ax2[0, 2].set_title("Mean Intensity Periphery/Cytoplasm Ratio")

sns.violinplot(
    data=df,
    x="RNAType",
    y="MeanRNAInPeriphery_GranuleCenters_Ratio",
    hue="Condition_RNA",
    ax=ax2[0, 3],
    inner="quart",
    palette=color_palette,
)
sns.swarmplot(
    data=df,
    x="RNAType",
    y="MeanRNAInPeriphery_GranuleCenters_Ratio",
    hue="Condition_RNA",
    ax=ax2[0, 3],
    color="k",
    alpha=0.5,
    dodge=True,
)
ax2[0, 3].set_title("Mean Intensity Periphery/Granule Centers Ratio")

sns.violinplot(
    data=df,
    x="RNAType",
    y="RNAGranuleOverlap_Fraction",
    hue="Condition_RNA",
    ax=ax2[1, 0],
    inner="quart",
    palette=color_palette,
)
sns.swarmplot(
    data=df,
    x="RNAType",
    y="RNAGranuleOverlap_Fraction",
    hue="Condition_RNA",
    ax=ax2[1, 0],
    color="k",
    alpha=0.5,
    dodge=True,
)
ax2[1, 0].set_title("RNA-Granule Overlap Fraction (Granule Area/ Total RNA Area)")

sns.violinplot(
    data=df,
    x="RNAType",
    y="RNAPunctaCount",
    hue="Condition_RNA",
    ax=ax2[1, 1],
    inner="quart",
    palette=color_palette,
)
sns.swarmplot(
    data=df,
    x="RNAType",
    y="RNAPunctaCount",
    hue="Condition_RNA",
    ax=ax2[1, 1],
    color="k",
    alpha=0.5,
    dodge=True,
)
ax2[1, 1].set_title("RNA Puncta Count")

sns.violinplot(
    data=df_per_sample,
    x="Condition",
    y="MeanGranuleSize",
    ax=ax2[1, 2],
    inner="quart",
    palette="Set2",
)
sns.swarmplot(
    data=df_per_sample,
    x="Condition",
    y="MeanGranuleSize",
    ax=ax2[1, 2],
    color="k",
    alpha=0.5,
    dodge=True,
)
ax2[1, 2].set_title("Mean Granule Size")


sns.violinplot(
    data=df_per_sample,
    x="Condition",
    y="GranuleCount",
    ax=ax2[1, 3],
    inner="quart",
    palette="Set2",
)
sns.swarmplot(
    data=df_per_sample,
    x="Condition",
    y="GranuleCount",
    ax=ax2[1, 3],
    color="k",
    alpha=0.5,
    dodge=True,
)
ax2[1, 3].set_title("Granule Count")

sns.despine(fig2)
plt.tight_layout()

fig2.savefig("output/quantifiaction.pdf")

fig_paper, ax_paper = plt.subplots(1, 1, figsize=(3, 4))
sns.violinplot(
    data=df,
    x="RNAType",
    y="MeanRNA_Cytoplasm_Granules_Ratio",
    hue="Condition",
    ax=ax_paper,
    inner="quart",
    split=True,
    # palette=color_palette,
)
sns.swarmplot(
    data=df,
    x="RNAType",
    y="MeanRNA_Cytoplasm_Granules_Ratio",
    hue="Condition",
    ax=ax_paper,
    color="k",
    alpha=0.85,
    dodge=True,
    size=3,
)
ax_paper.set_title("Mean Intensity (whole) Granules/Cytoplasm Ratio")
ax_paper.set_ylabel(r"$\frac{\overline{I}_{Cytoplasm}}{\overline{I}_{Granule}}$")
plt.tight_layout()
sns.despine(fig_paper)
fig_paper.savefig("output/quantification_paper.pdf")

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 6,  # default text
        "axes.titlesize": 8,  # subplot titles
        "axes.labelsize": 6,  # x/y axis labels
        "xtick.labelsize": 6,  # x tick labels
        "ytick.labelsize": 6,  # y tick labels
        "legend.fontsize": 6,  # legend
        "figure.titlesize": 8,  # suptitle
    }
)

color_palette2 = {
    "DMSO Control": "#888888",
    "CHX Treated": "#224644",
    "PatA": "#566038",
}
fig_paper_newstyle, ax_paper_newstyle = plt.subplots(1, 1, figsize=(1.63, 3))

sns.stripplot(
    data=df,
    x="RNAType",
    y="MeanRNA_Cytoplasm_Granules_Ratio",
    hue="Condition",
    dodge=True,
    palette=color_palette2,
    alpha=0.4,
    size=3,
    ax=ax_paper_newstyle,
    jitter=0.3,
)
sns.pointplot(
    data=df,
    x="RNAType",
    y="MeanRNA_Cytoplasm_Granules_Ratio",
    hue="Condition",
    dodge=0.4,
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
    ax=ax_paper_newstyle,
)
ax_paper_newstyle.set_title("mRNA relocalization - CHX", fontweight="bold")
ax_paper_newstyle.set_ylabel(
    r"Mean Intensity Ratio (Cytoplasm/Granule)", fontweight="bold"
)
ax_paper_newstyle.set_xticklabels(
    ax_paper_newstyle.get_xticklabels(), fontstyle="italic"
)
ax_paper_newstyle.set_xlabel("")
plt.tight_layout()
sns.despine()
fig_paper_newstyle.savefig(
    "output/quantification_paper_newstyle.pdf", dpi=600, transparent=True
)

for feature in [
    "MeanRNAInGranules_Cytoplasm_Ratio",
    "MeanRNAInGranuleCenters_Cytoplasm_Ratio",
    "MeanRNAInPeriphery_Cytoplasm_Ratio",
    "MeanRNAInPeriphery_GranuleCenters_Ratio",
    "RNAGranuleOverlap_Fraction",
    "RNAPunctaCount",
]:
    for rna_type in ["nanos", "tdrd7"]:
        significance_testing.append(get_signifcances(df, feature, rna_type))

df_significance = pd.DataFrame(significance_testing)
df_significance.to_csv("output/statistical_testing_results.csv", index=False)


def get_representative_ids(df, feature, rep_min=40, rep_max=60):
    """
    For a given feature returns the IDs of samples that are representative all conditions, meaning they are around the median of the distribution for that feature. This can be used to select samples for visualizations that are representative of the overall trends in the data.
    params:
    df: DataFrame containing results of an expeiment
    feature: Feature, meaning column, of the dataframe in the
    rep_min: lower percentile threshold for selecting representative samples
    rep_max: upper percentile threshold for selecting representative samples
    """
    condition_dict = {}
    for condition in df["Condition"].unique():
        representative_ids = []
        for rna_type in df["RNAType"].unique():
            sub_df = df.query("RNAType == @rna_type and Condition == @condition")
            min_val = np.percentile(sub_df[feature].values, rep_min)
            max_val = np.percentile(sub_df[feature].values, rep_max)
            if len(representative_ids) == 0:
                # print(sub_df.head(10))
                # print(sub_df.columns)
                representative_ids = sub_df.query(f"@min_val <= {feature} <= @max_val")[
                    "UID"
                ].values
            else:
                representative_ids = np.intersect1d(
                    representative_ids,
                    sub_df.query(f"@min_val <= {feature} <= @max_val")["UID"].values,
                )

            condition_dict[condition] = representative_ids

    return condition_dict


with open("output/representative_samples.txt", "w") as f:
    for feature in [
        "MeanRNAInGranules_Cytoplasm_Ratio",
        "MeanRNAInGranuleCenters_Cytoplasm_Ratio",
        "MeanRNAInPeriphery_Cytoplasm_Ratio",
        "MeanRNAInPeriphery_GranuleCenters_Ratio",
        "RNAGranuleOverlap_Fraction",
        "RNAPunctaCount",
    ]:
        f.write(f"Feature: {feature}\n")
        representative_ids = get_representative_ids(df, feature)
        for condition, ids in representative_ids.items():
            f.write(f"  Condition: {condition}\n")
            for uid in ids:
                f.write(f"    Representative Sample ID: {uid}\n")


df.to_csv("output/quantification_results.csv", index=False)
