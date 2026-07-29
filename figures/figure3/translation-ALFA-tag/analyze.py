import pandas as pd
from segmenter import SegmentationModel
from pathlib import Path
import matplotlib.pyplot as plt
from skimage.filters import gaussian
from skimage.measure import label
import scipy.ndimage as ndi
import numpy as np
import tifffile as tiff
from utils import (
    select_center_object,
    lsm_pixel_size,
    get_membrane_mask,
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

fig_contact_sheet, ax_contact_sheet = plt.subplots(
    len(list_of_samples),
    6,
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

    # Write masks to file
    if not recompute_masks or not (p / "masks" / "cell.tif").exists():
        mask_dir = p / "masks"
        mask_dir.mkdir(exist_ok=True, parents=True)

        tiff.imwrite(mask_dir / "cell.tif", cell_seg.astype(np.uint8))
        tiff.imwrite(mask_dir / "nucleus.tif", nucleus_seg.astype(np.uint8))
        tiff.imwrite(mask_dir / "granule.tif", granule_seg.astype(np.uint8))

    plot_image(img_vasa, ax[0], colormap="cyan")
    ax[0].contour(granule_seg > 0, colors="white", linewidths=0.5)
    plot_image(img_mem, ax[1])
    ax[1].contour(cell_seg > 0, colors="white", linewidths=0.5)
    plot_image(
        img_mCh, ax[2], colormap="magenta", vmin=0, vmax=np.percentile(img_mCh, 99.999)
    )
    ax[2].contour(nucleus_seg > 0, colors="white", linewidths=0.5)
    plot_image(img_alfatag, ax[3], vmin=0, vmax=np.percentile(img_alfatag, 87))
    plot_merge(
        {
            "cyan": img_vasa,
            "magenta": img_mCh,
        },
        ax[4],
        vmins={"cyan": 0, "magenta": 0},
        vmaxs={
            "cyan": np.percentile(img_vasa, 99.999),
            "magenta": np.percentile(img_mCh, 80),
        },
    )
    for a in ax:
        a.axis("off")


fig_contact_sheet.savefig(output_dir / "contact_sheet.pdf", bbox_inches="tight")
