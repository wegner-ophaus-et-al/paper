import tifffile as tiff
import matplotlib.pyplot as plt
from skimage import measure


def delete_outside_segmentations(segmentation_of_interest, confining_segmentation):
    out_seg = segmentation_of_interest.copy()
    out_seg[~confining_segmentation.astype(bool)] = 0
    return out_seg


def process_sample(
    sample_folder,
    segmenter_instance,
    resegmentation=False,
    ax=None,
):

    out_data_dict = {}

    uid, name = sample_folder.name.split("__", 1)
    construct_id = sample_folder.parent.name
    condition = (
        "control"
        if "control" in name.lower() or "ctrl" in name.lower()
        else "heatshock"
    )

    out_data_dict.update(
        {
            "uid": uid,
            "sample_name": name,
            "construct_id": construct_id,
            "condition": condition,
        }
    )

    int_img = {
        "granulito": tiff.imread(sample_folder / "images" / "gra.tif"),
        "nls": tiff.imread(sample_folder / "images" / "nls.tif"),
        "stress_marker": tiff.imread(sample_folder / "images" / "stress.tif"),
    }

    lbl_img = {
        "cell": measure.label(tiff.imread(sample_folder / "masks" / "cell.tif")),
    }

    if not (sample_folder / "masks" / "granules.tif").exists() or resegmentation:
        gra_segmentation = segmenter_instance.segment(int_img["granulito"])
        tiff.imwrite(
            sample_folder / "masks" / "granules.tif", gra_segmentation.astype("uint8")
        )

    else:
        gra_segmentation = tiff.imread(sample_folder / "masks" / "granules.tif")

    lbl_img["granules"] = delete_outside_segmentations(
        gra_segmentation, lbl_img["cell"]
    )
    lbl_img["cytoplasm"] = lbl_img["cell"].astype(bool) & ~lbl_img["granules"].astype(
        bool
    )

    out_data_dict.update(
        {
            "total_cell_area": (lbl_img["cell"] > 0).sum(),
            "total_granules_area": (lbl_img["granules"] > 0).sum(),
            "total_cytoplasm_area": (lbl_img["cytoplasm"] > 0).sum(),
            "mean_int_granulito_cell": int_img["granulito"][lbl_img["cell"] > 0].mean(),
            "mean_int_granulito_granules": int_img["granulito"][
                lbl_img["granules"] > 0
            ].mean(),
            "mean_int_granulito_cytoplasm": int_img["granulito"][
                lbl_img["cytoplasm"]
            ].mean(),
            "mean_int_stress_cell": int_img["stress_marker"][
                lbl_img["cell"] > 0
            ].mean(),
            "mean_int_stress_granules": int_img["stress_marker"][
                lbl_img["granules"] > 0
            ].mean(),
            "mean_int_stress_cytoplasm": int_img["stress_marker"][
                lbl_img["cytoplasm"] > 0
            ].mean(),
        }
    )

    if ax is None:
        fig, ax = plt.subplots(1, 3, figsize=(10, 3))
        print("You're in the wrong hole")
    else:
        fig = None

    ax[0].imshow(int_img["granulito"], cmap="gray")
    ax[0].contour(lbl_img["granules"], colors="r", linewidths=0.5)
    ax[0].set_title("Granulito")
    ax[0].text(
        -0.05,
        0.5,  # x slightly outside the left edge, y centered
        f"{uid}-{condition}",
        transform=ax[0].transAxes,
        rotation=90,
        va="center",
        ha="right",
    )

    ax[1].imshow(int_img["stress_marker"], cmap="gray")
    ax[1].contour(lbl_img["granules"], colors="r", linewidths=0.5)
    ax[1].set_title("Stress Marker")

    ax[2].imshow(int_img["nls"], cmap="gray")
    ax[2].contour(lbl_img["cell"], colors="r", linewidths=0.5)
    ax[2].set_title("NLS")

    if fig is not None:
        plt.savefig(sample_folder / "visualization.png", dpi=300)
        plt.close(fig)

    return out_data_dict
