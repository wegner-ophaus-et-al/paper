import pandas as pd
from segmenter import SegmentationModel
from pathlib import Path
import matplotlib.pyplot as plt
from skimage.filters import gaussian
import scipy.ndimage as ndi
import numpy as np
import tifffile as tiff
from utils import select_center_object, lsm_pixel_size
from filtering import segment_spots, dog_filter
from plotlib import plot_merge, plot_image, scale_bar
import seaborn as sns


recompute_masks = False


root = Path("/Volumes/Kur/paper_data_sorted/tmd-nanos/2026-07-27_TMD-nanos_confo_6hpf")
output_dir = Path(__file__).parent / "output"
output_dir.mkdir(exist_ok=True, parents=True)

list_of_samples = [
    p for p in root.iterdir() if p.is_dir() and not p.name.startswith(".")
]

sm = SegmentationModel(path=None, model_type="vit_l_lm", upsampling_factor=1)


contact_sheet_size = 4

fig_contact_sheet, ax_contact_sheet = plt.subplots(
    len(list_of_samples),
    6,
    figsize=(6 * contact_sheet_size, len(list_of_samples) * contact_sheet_size * 0.9),
)

# Set Dnd1 spot segmentation params
sig = 2.5
k_value = 6


results = []
for ax, p in zip(ax_contact_sheet, list_of_samples):
    results_dict = {}

    uid, sample_name = p.name.split("__")
    condition_string, cell_number = sample_name.split("-")

    injected_components = [c.upper() for c in condition_string.split("_")]

    condition = ""
    if "A709" in injected_components and "F125" in injected_components:
        condition = "wt"
    elif "F125" in injected_components and "F649" in injected_components:
        condition = "tmd-n_stuffer"
    elif "B008" in injected_components and "F649" in injected_components:
        condition = "tmd-n_d"
    elif "B008" in injected_components and "D378" in injected_components:
        condition = "n_d"
    else:
        raise KeyError(f"{condition_string} could not be assinged to any condition")

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

    # Dnd channel
    img_dnd = img[1]
    img_dnd_filteres = dog_filter(img_dnd, sigma=sig)
    dnd_spots, spot_number, threshold_value = segment_spots(
        img[1], sigma=sig, k=k_value
    )
    spots_area = np.sum(dnd_spots > 0) * pixel_size**2

    # Membrane channel
    img_mem = img[2].copy()
    if not recompute_masks and (p / "masks" / "cell.tif").exists():
        cell_seg = tiff.imread(p / "masks" / "cell.tif")

    else:
        cell_seg = sm.segment(gaussian(img_mem, sigma=1.5))

    cell_seg = select_center_object(cell_seg)
    img_distance_from_membrane = ndi.distance_transform_edt(cell_seg > 0)
    img_distance_from_membrane = img_distance_from_membrane * pixel_size

    # Get spot distances
    img_dfm_spots = img_distance_from_membrane.copy()
    img_dfm_spots[dnd_spots < 1] = 0
    spot_distances = img_dfm_spots[dnd_spots > 0].flatten()

    results_dict.update(
        {
            "vasa_mean_intensity": np.mean(img_vasa),
            "membrane_mean_intensity": np.mean(img_mem),
            "cell_area": np.sum(cell_seg > 0) * pixel_size**2,
            "dnd_mean_intensity": np.mean(img_dnd),
            "spot_number": spot_number,
            "spots_area": spots_area,
            "spots_threshold": threshold_value,
            "spot_distances_mean": np.mean(spot_distances),
            "spot_distances": spot_distances,
        }
    )

    plot_merge(
        {
            "cyan": img_vasa,
            "magenta": img_mem,
            "yellow": img_dnd,
        },
        ax=ax[0],
    )
    plot_image(img_vasa, ax[1], colormap="cyan")
    ax[1].contour(cell_seg, colors="white", linewidths=0.5)
    plot_image(img_mem, ax[2], colormap="magenta")
    plot_image(img_dnd, ax[3], colormap="yellow")
    plot_image(img_dnd_filteres, ax[4], colormap="yellow")
    ax[4].contour(dnd_spots, colors="white", linewidths=0.5)
    ax[4].imshow(img_dnd_filteres, cmap="magma")
    ax[4].set_title(f"{uid}-{condition}")
    ax[5].imshow(img_dfm_spots, cmap="magma")

    for sb_ax_idx in range(5):
        scale_bar(ax[sb_ax_idx], length=5, thickness=1, pixel_size=pixel_size)

    # Write masks to file
    mask_dir = p / "masks"
    mask_dir.mkdir(exist_ok=True, parents=True)

    tiff.imwrite(mask_dir / "cell.tif", cell_seg.astype(np.uint8))
    tiff.imwrite(mask_dir / "spots.tif", dnd_spots.astype(np.uint8))

    for a in ax.flatten():
        a.axis("off")

    results.append(results_dict)


fig_contact_sheet.tight_layout()
fig_contact_sheet.savefig(output_dir / "contact_sheet.pdf")
plt.close(fig_contact_sheet)

df_cell = pd.DataFrame(results)
distances_by_condition = []

for res in results:
    for dist in res["spot_distances"]:
        distances_by_condition.append({"condition": res["condition"], "distance": dist})
df_dist_by_cond = pd.DataFrame(distances_by_condition)

sns.set_style("ticks")
plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 6,  # default text
        "axes.titlesize": 8,  # subplot titles
        "axes.labelsize": 6,  # x/y axis labels
        "axes.labelweight": "bold",  # x/y axis labels
        "xtick.labelsize": 6,  # x tick labels
        "ytick.labelsize": 6,  # y tick labels
        "legend.fontsize": 6,  # legend
        "figure.titlesize": 8,  # suptitle
    }
)

fig_kde, ax_kde = plt.subplots(figsize=(5, 3))
sns.kdeplot(
    data=df_dist_by_cond, x="distance", hue="condition", ax=ax_kde, common_norm=False
)
ax_kde.set_xlabel("Distance from membrane (µm)")
ax_kde.set_title("Dnd1 spot distances from membrane")

fig_kde.tight_layout()
sns.despine()
fig_kde.savefig(output_dir / "spot_distances_kde.pdf")
plt.close(fig_kde)

ncols = 8
nrows = 1
width = 1.4
height = 2.5

fig_cell_features, axs = plt.subplots(
    nrows=nrows, ncols=ncols, figsize=(ncols * width, nrows * height)
)

for idax, feature in enumerate(
    [
        "vasa_mean_intensity",
        "membrane_mean_intensity",
        "cell_area",
        "dnd_mean_intensity",
        "spot_number",
        "spots_area",
        "spots_threshold",
        "spot_distances_mean",
    ]
):
    sns.stripplot(
        data=df_cell,
        y=feature,
        hue="condition",
        # hue_order=["full_mix", "no_tdrd7a"],
        # palette=color_palette,
        dodge=0.4,
        alpha=0.4,
        size=2,
        ax=axs[idax],
        jitter=0.2,
    )
    # sns.pointplot(
    #     data=df_cell,
    #     y=feature,
    #     hue="condition",
    #     # hue_order=["full_mix", "no_tdrd7a"],
    #     dodge=0.4,
    #     ax=axs[idax],
    #     errorbar="sd",  # standard error
    #     estimator="median",  # or "mean"
    #     capsize=0.075,
    #     linestyle="none",
    #     markersize=10,
    #     marker="_",
    #     err_kws=dict(linewidth=0.5, color="black"),
    #     markeredgewidth=1,
    #     palette="dark:black",
    #     # zorder=5,
    # )

    axs[idax].set_title(feature.replace("_", " ").title(), fontweight="bold")
    axs[idax].xaxis.label.set_fontweight("bold")

fig_cell_features.tight_layout()
sns.despine()
fig_cell_features.savefig(output_dir / "cell_features.pdf")
plt.close(fig_cell_features)
