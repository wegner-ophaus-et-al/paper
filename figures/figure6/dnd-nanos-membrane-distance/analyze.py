from segmenter import SegmentationModel
from pathlib import Path
import matplotlib.pyplot as plt
from skimage.filters import gaussian
import scipy.ndimage as ndi
import numpy as np
import tifffile as tiff
from utils import select_center_object, lsm_pixel_size
from filtering import segment_spots, dog_filter
from plotlib import plot_merge, plot_image

recompute_masks = False


root = Path("/Volumes/Kur/paper_data_sorted/tmd-nanos/2026-07-27_TMD-nanos_confo_6hpf")
output_dir = Path(__file__).parent / "output"
output_dir.mkdir(exist_ok=True, parents=True)

list_of_samples = [
    p for p in root.iterdir() if p.is_dir() and not p.name.startswith(".")
]

sm = SegmentationModel(path=None, model_type="vit_l_lm", upsampling_factor=1)


contact_sheet_size = 2

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
        condition = "n-d"
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
    spots_area = np.sum(dnd_spots > 0)

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
            "dnd_mean_intensity": np.mean(img_dnd),
            "membrane_mean_intensity": np.mean(img_mem),
            "spot_number": spot_number,
            "spots_area": spots_area,
            "spots_threshold": threshold_value,
            "spot_distances_mean": np.mean(spot_distances),
            "spot_distances": spot_distances,
        }
    )

    # plot_merge(
    #     {
    #         "cyan": img_vasa,
    #         "magenta": img_mem,
    #         "yellow": img_dnd_filteres,
    #     },
    #     ax=ax,
    # )
    plot_image(img_vasa, ax[1], colormap="cyan")
    ax[1].contour(cell_seg, colors="white", linewidths=0.5)
    plot_image(img_mem, ax[2], colormap="magenta")
    plot_image(img_dnd, ax[3], colormap="yellow")
    plot_image(img_dnd_filteres, ax[4], colormap="yellow")
    ax[4].contour(dnd_spots, colors="white", linewidths=0.5)
    ax[4].imshow(img_dnd_filteres, cmap="magma")
    ax[4].set_title(f"{uid}-{condition}")
    ax[5].imshow(img_dfm_spots, cmap="magma")

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
