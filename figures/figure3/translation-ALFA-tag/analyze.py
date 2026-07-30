import pandas as pd
from segmenter import SegmentationModel
from pathlib import Path
import matplotlib.pyplot as plt
from skimage.filters import gaussian
from skimage.measure import regionprops, label
import scipy.ndimage as ndi
import numpy as np
import tifffile as tiff
from utils import (
    select_center_object,
    lsm_pixel_size,
    get_ring_mask,
    confine_mask_to_cell,
    statistical_analysis,
    get_stars,
)
from filtering import segment_spots, dog_filter
from plotlib import plot_merge, plot_image, scale_bar
from pathlib import Path

# from representative_images import export_representative_uids
import seaborn as sns


root = Path("/Volumes/Kur/paper_data_sorted/ALFAtag/19xALFAtag_confo_KT/all_data")
# root = Path("/Users/icb_remote/Documents/JW/py/data/19xALFAtag_confo_KT/all_data/")
output_dir = Path(__file__).parent / "output"
output_dir.mkdir(exist_ok=True, parents=True)

list_of_samples = [
    p for p in root.iterdir() if p.is_dir() and not p.name.startswith(".")
]


recompute_masks = False

# sm = SegmentationModel(path=None, model_type="vit_l_lm", upsampling_factor=1)
sm = None

sm_g = None
# sm_g = SegmentationModel(
#     path=Path(
#         "/Users/icb_remote/Documents/JW/py/packages/gamgee/models/msam/sam_granules_refined_up3_34917658/"
#     ),
#     mpdel_type="vit_l_lm",
#     upsampling_factor=3,
# )

contact_sheet_size = 4

###############
sigma = 2
k_value = 6
safety_barrier = 1
###############


fig_contact_sheet, ax_contact_sheet = plt.subplots(
    len(list_of_samples),
    7,
    figsize=(6 * contact_sheet_size, len(list_of_samples) * contact_sheet_size * 0.9),
)

results = []
for ax, p in zip(ax_contact_sheet, list_of_samples):
    results_dict = {}

    uid, sample_name = p.name.split("__")

    condition = ""
    if "pata" in sample_name.lower():
        condition = "PatA"
    else:
        condition = "Control"

    lsm_dir = p / "original_raw"
    lsm_path = next(lsm_dir.glob("*.lsm"), None)
    pixel_size = lsm_pixel_size(lsm_path)

    results_dict.update(
        {
            "uid": uid,
            "name": sample_name,
            "condition": condition,
            "lsm_file": lsm_path.name,
            "pixel_size": pixel_size,
        }
    )

    img = tiff.imread(lsm_path)

    # Vasa channel
    img_vasa = img[0]
    if not recompute_masks and (p / "masks" / "granule.tif").exists():
        granule_seg = tiff.imread(p / "masks" / "granule.tif")
    else:
        granule_seg = sm_g.segment(gaussian(img_vasa, sigma=1.5))

    # Membrane channel
    img_mem = img[1].copy()
    if not recompute_masks and (p / "masks" / "cell.tif").exists():
        cell_seg = tiff.imread(p / "masks" / "cell.tif")
    else:
        cell_seg = sm.segment(gaussian(img_mem, sigma=1.5))

    cell_seg = select_center_object(cell_seg)
    # Target expression channel
    img_mCh = img[2]

    # ALFAtag
    img_alfatag = img[3]
    if not recompute_masks and (p / "masks" / "nucleus.tif").exists():
        nucleus_seg = tiff.imread(p / "masks" / "nucleus.tif")
    else:
        nucleus_seg = sm.segment(gaussian(img_alfatag, sigma=1.5))

    # Delete outside cell objects
    nucleus_seg = confine_mask_to_cell(nucleus_seg, cell_seg)
    granule_seg = confine_mask_to_cell(granule_seg, cell_seg)

    # Generate masks
    granule_periphery_mask = get_ring_mask(
        granule_seg, thickness=safety_barrier, pixel_size=pixel_size
    )
    nucleus_periphery_mask = get_ring_mask(
        nucleus_seg, thickness=safety_barrier, pixel_size=pixel_size
    )
    cytoplasm_mask = (
        cell_seg
        - nucleus_seg
        - granule_seg
        - granule_periphery_mask
        - nucleus_periphery_mask
    )

    # Generate ALFA-tag spots
    alfa_spots, num_spots, _ = segment_spots(
        img_alfatag, sigma=sigma, k=k_value, min_area=4, offset=0.0
    )

    # Remove spots outside the cell and inside the nucleus+safety barrier
    alfa_spots = confine_mask_to_cell(alfa_spots, cell_seg)
    alfa_spots = confine_mask_to_cell(alfa_spots, nucleus_seg == 0)
    alfa_spots = confine_mask_to_cell(alfa_spots, nucleus_periphery_mask == 0)

    # Write masks to file
    if not recompute_masks or not (p / "masks" / "cell.tif").exists():
        mask_dir = p / "masks"
        mask_dir.mkdir(exist_ok=True, parents=True)

        tiff.imwrite(mask_dir / "cell.tif", cell_seg.astype(np.uint8))
        tiff.imwrite(mask_dir / "nucleus.tif", nucleus_seg.astype(np.uint8))
        tiff.imwrite(mask_dir / "granule.tif", granule_seg.astype(np.uint8))

    # Plotting images to contanct sheet
    plot_image(img_vasa, ax[0], colormap="cyan")
    ax[0].contour(granule_seg > 0, colors="white", linewidths=0.5)
    plot_image(img_mem, ax[1])
    ax[1].contour(cell_seg > 0, colors="white", linewidths=0.5)
    plot_image(
        img_mCh, ax[2], colormap="yellow", vmin=0, vmax=np.percentile(img_mCh, 99.999)
    )
    ax[2].contour(nucleus_seg > 0, colors="white", linewidths=0.5)
    plot_image(img_alfatag, ax[3], vmin=0, vmax=np.percentile(img_alfatag, 95))
    plot_merge(
        {
            "cyan": img_vasa,
            "magenta": img_mCh,
        },
        ax[4],
        vmins={"cyan": 0, "magenta": 0},
        vmaxs={
            "cyan": np.percentile(img_vasa, 99.999),
            "magenta": np.percentile(img_alfatag, 95),
        },
    )
    ax[5].imshow(np.zeros(cell_seg.shape), cmap="gray")
    ax[5].contour(cytoplasm_mask > 0, colors="#FF00FF", linewidths=0.5)
    ax[5].contour(nucleus_seg > 0, colors="#00FFFF", linewidths=0.5)
    ax[5].contour(granule_seg > 0, colors="#ffff00", linewidths=0.5)
    ax[6].imshow(alfa_spots.astype("float32") * (255 / np.max(alfa_spots)), cmap="gray")
    for a in ax:
        a.axis("off")

    # Measurement and statistical analysis
    nucleus_target_expression = np.mean(img_mCh[nucleus_seg > 0])

    mean_alfa_cytoplasm = np.mean(img_alfatag[cytoplasm_mask > 0])
    sum_alfa_cytoplasm = np.sum(img_alfatag[cytoplasm_mask > 0])
    cv_alfa_cytoplasm = np.std(img_alfatag[cytoplasm_mask > 0]) / mean_alfa_cytoplasm

    mean_alfa_granules = np.mean(img_alfatag[granule_seg > 0])
    sum_alfa_granules = np.sum(img_alfatag[granule_seg > 0])
    cv_alfa_granules = np.std(img_alfatag[granule_seg > 0]) / mean_alfa_granules

    mean_alfa_granule_periphery = np.mean(img_alfatag[granule_periphery_mask > 0])
    sum_alfa_granule_periphery = np.sum(img_alfatag[granule_periphery_mask > 0])
    cv_alfa_granule_periphery = (
        np.std(img_alfatag[granule_periphery_mask > 0]) / mean_alfa_granule_periphery
    )

    alfa_spots_labeled = label(alfa_spots)
    number_of_spots = np.max(alfa_spots_labeled)
    area_of_spots = np.sum(alfa_spots > 0)
    sum_alfa_spots = np.sum(img_alfatag[alfa_spots > 0])

    granule_distance_map = ndi.distance_transform_edt(granule_seg == 0) * pixel_size
    distances_to_granules = []
    for region in regionprops(alfa_spots_labeled):
        distances_to_granules.append(
            granule_distance_map[tuple(np.round(region.centroid).astype(int))]
        )
    mean_distance_to_granules = (
        np.mean(distances_to_granules) if distances_to_granules else np.nan
    )

    results_dict.update(
        {
            "nucleus_target_expression": nucleus_target_expression,
            "mean_alfa_cytoplasm": mean_alfa_cytoplasm,
            "sum_alfa_cytoplasm": sum_alfa_cytoplasm,
            "cv_alfa_cytoplasm": cv_alfa_cytoplasm,
            "mean_alfa_granules": mean_alfa_granules,
            "sum_alfa_granules": sum_alfa_granules,
            "cv_alfa_granules": cv_alfa_granules,
            "mean_alfa_granule_periphery": mean_alfa_granule_periphery,
            "sum_alfa_granule_periphery": sum_alfa_granule_periphery,
            "cv_alfa_granule_periphery": cv_alfa_granule_periphery,
            "number_of_spots": number_of_spots,
            "area_of_spots": area_of_spots,
            "sum_intensity_alfa_spots": sum_alfa_spots,
            "mean_distance_to_granules": mean_distance_to_granules,
        }
    )

    results.append(results_dict)


fig_contact_sheet.savefig(output_dir / "contact_sheet.pdf", bbox_inches="tight")


df = pd.DataFrame(results)
df["MeanRNA_Cytoplasm_Granules_Ratio"] = (
    df["mean_alfa_cytoplasm"] / df["mean_alfa_granules"]
)
df["SumRNA_Cytoplasm_Granules_Ratio"] = (
    df["sum_alfa_cytoplasm"] / df["sum_alfa_granules"]
)


# Statistical analysis: Control vs PatA for every measured feature
all_features = [
    "nucleus_target_expression",
    "mean_alfa_cytoplasm",
    "sum_alfa_cytoplasm",
    "cv_alfa_cytoplasm",
    "mean_alfa_granules",
    "sum_alfa_granules",
    "cv_alfa_granules",
    "mean_alfa_granule_periphery",
    "sum_alfa_granule_periphery",
    "cv_alfa_granule_periphery",
    "number_of_spots",
    "area_of_spots",
    "sum_intensity_alfa_spots",
    "mean_distance_to_granules",
    "MeanRNA_Cytoplasm_Granules_Ratio",
    "SumRNA_Cytoplasm_Granules_Ratio",
]

stats_results = []
for feature in all_features:
    control_vals = df.loc[df["condition"] == "Control", feature].dropna().to_numpy()
    pata_vals = df.loc[df["condition"] == "PatA", feature].dropna().to_numpy()
    test_type, statistic, p_value = statistical_analysis(control_vals, pata_vals)
    stars = get_stars(p_value)
    stats_results.append(
        {
            "feature": feature,
            "test_type": test_type,
            "n_control": len(control_vals),
            "n_pata": len(pata_vals),
            "statistic": statistic,
            "p_value": p_value,
            "stars": stars,
        }
    )

stats_df = pd.DataFrame(stats_results)
stats_df.to_csv(output_dir / "stats_results.csv", index=False)
stars_by_feature = dict(zip(stats_df["feature"], stats_df["stars"]))


sns.set_style("ticks")

plt.rcParams.update(
    {
        "font.size": 6,  # default text
        "font.family": "Arial",
        "axes.titlesize": 8,  # subplot titles
        "axes.labelsize": 6,  # x/y axis labels
        "axes.labelweight": "bold",  # x/y axis labels
        "xtick.labelsize": 6,  # x tick labels
        "ytick.labelsize": 6,  # y tick labels
        "legend.fontsize": 6,  # legend
        "figure.titlesize": 8,  # suptitle
    }
)

color_palette = {
    "Control": "#888888",
    "CHX": "#224644",
    "PatA": "#566038",
}

ncols_basic = 3
nrows_basic = 3
width = 1.4
height = 2.5

fig_basic, ax_basic = plt.subplots(
    nrows=nrows_basic,
    ncols=ncols_basic,
    figsize=(ncols_basic * width, nrows_basic * height),
)
for feature, ax in zip(
    [
        "mean_alfa_cytoplasm",
        "sum_alfa_cytoplasm",
        "cv_alfa_cytoplasm",
        "mean_alfa_granules",
        "sum_alfa_granules",
        "cv_alfa_granules",
        "mean_alfa_granule_periphery",
        "sum_alfa_granule_periphery",
        "cv_alfa_granule_periphery",
    ],
    ax_basic.flatten(),
):
    sns.stripplot(
        data=df,
        # x="RNAType",
        y=feature,
        hue="condition",
        dodge=True,
        palette=color_palette,
        alpha=0.4,
        size=3,
        ax=ax,
        jitter=0.3,
    )
    sns.pointplot(
        data=df,
        # x="RNAType",
        y=feature,
        hue="condition",
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
        ax=ax,
    )
    stars = stars_by_feature.get(feature, "")
    ax.set_title(
        f"{str(feature).replace('_', ' ').title()} ({stars})", fontweight="bold"
    )
    ax.set_ylabel(str(feature).replace("_", " "), fontweight="bold")
fig_basic.tight_layout()
sns.despine()
fig_basic.savefig(output_dir / "basic_features.pdf", bbox_inches="tight")


ncols_spots = 3
nrows_spots = 2

fig_spots, ax_spots = plt.subplots(
    nrows=nrows_spots,
    ncols=ncols_spots,
    figsize=(ncols_spots * width, nrows_spots * height),
)
for ax_2, feature in zip(
    ax_spots.flatten(),
    [
        "number_of_spots",
        "area_of_spots",
        "sum_intensity_alfa_spots",
        "mean_distance_to_granules",
        "nucleus_target_expression",
        "SumRNA_Cytoplasm_Granules_Ratio",
    ],
):
    sns.stripplot(
        data=df,
        y=feature,
        hue="condition",
        dodge=True,
        palette=color_palette,
        alpha=0.4,
        size=3,
        ax=ax_2,
        jitter=0.3,
    )
    sns.pointplot(
        data=df,
        y=feature,
        hue="condition",
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
        ax=ax_2,
    )
    stars = stars_by_feature.get(feature, "")
    ax_2.set_title(
        f"{str(feature).replace('_', ' ').title()} ({stars})", fontweight="bold"
    )
    ax_2.set_ylabel(str(feature).replace("_", " "), fontweight="bold")
fig_spots.tight_layout()
sns.despine()
fig_spots.savefig(output_dir / "spots_features.pdf", bbox_inches="tight")
